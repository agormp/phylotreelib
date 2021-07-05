"""Classes and methods for analyzing, manipulating, and building phylogenetic trees"""
# Anders Gorm Pedersen
# Section for Bioinformatics, DTU Health Technology, Technical University of Denmark
# agpe@dtu.dk
# Requires Python 3.0+

import copy
import functools
import itertools
import math
import random
import re
import statistics
import sys
from io import StringIO
import numpy as np

#############################################################################################
#############################################################################################
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
##                  check that this behaves as expected wrt performance (and remove Globals perhaps)
##
##        (2) Although there is some extra overhead in using classes to emulate structs,
##            using dicts instead does not make a big difference performancewide (I tried).
##
##        (3) Treestrings could be parsed recursively in fewer lines, but this works a lot
##            less efficiently than the iterated version currently in Tree.from_string()
##            (I tested it).
##
#############################################################################################
#############################################################################################

#############################################################################################
#############################################################################################
#
# Various functions used by methods, that do not fit neatly in any class

def escape_metachars(text, metachars=".^$*+?{}[]\|()"):
    """Add backslashes to escape metachars in input string."""

    newtxt = ""
    for char in text:
        if char in metachars:
            newtxt += "\\" + char
        else:
            newtxt += char
    return newtxt

###################################################################################################

def remove_comments(text, leftdelim, rightdelim=None):
    """Takes input string and strips away commented text, delimited by 'leftdelim' and 'rightdelim'.
        Also deals with nested comments."""

    # NOTE: only deals with block comments at present
    # Sanity checks: delimiters can not be identical, and one cannot be substring of the other
    if leftdelim == rightdelim:
        raise Exception("Left and right delimiters are identical")
    if leftdelim in rightdelim:
        raise Exception("Left delimiter is substring of right delimiter")
    if rightdelim in leftdelim:
        raise Exception("Right delimiter is substring of left delimiter")

    # Preprocess delims for use in re etc
    leftdelim = escape_metachars(leftdelim)
    rightdelim = escape_metachars(rightdelim)

    # Construct sorted list of tuples of the form [(0, 'start'), (5, 'stop'), (7, 'start'), ...]
    delimlist = [(delim.start(), "start") for delim in re.finditer(leftdelim, text)]
    # If text contains no starts (=> no comments): return un-altered text
    if not delimlist:
        return text
    delimlist.extend([(delim.start(), "stop") for delim in re.finditer(rightdelim, text)])
    delimlist.sort()

    # Resolve issues with overlapping delimiters:
    tmplist = [delimlist[0]]
    for i in range(1, len(delimlist)):
        start1 = delimlist[i - 1][0]
        start2 = delimlist[i][0]
        if start2 > start1 + len(leftdelim) - 1:    # If  next item does not overlap previous item
            tmplist.append(delimlist[i])
    delimlist = tmplist

    # Traverse text; along the way copy text not inside comment-delimiter pairs.
    # Use stack ("unmatched_starts") to keep track of nesting
    offset = len(rightdelim) - 1
    unmatched_starts = 0
    prevpos = 0
    processed_text = ""
    for (start, delim) in delimlist:
        if delim == "start":
            unmatched_starts += 1
            if unmatched_starts == 1:                   # Beginning of new comment region
                processed_text += text[prevpos:start]
        elif delim == "stop":
            unmatched_starts -= 1
            if unmatched_starts == 0:                   # End of comment region
                prevpos = start + offset
            elif unmatched_starts == -1:                # Error: more right delims than left delims
                raise Exception("Unmatched end-comment delimiter: {}".format(text[prevpos-5:prevpos+5]))

    # Add parts of original text to the right of rightdelim (if present)
    if prevpos < len(text):
        processed_text += text[prevpos:]
    return processed_text

#############################################################################################
#############################################################################################

class Globals():
    """Class containing globally used functions and labels."""

    # I'm not convinced this is the way to go. Module instead?"""
    # Global repository for bipartitions, to avoid redundant saving in topologies etc.
    biparts = {}

#############################################################################################
#############################################################################################

class Branchstruct():
    """Class that emulates a struct. Keeps branch-related info"""

    # Python note: perhaps replace with dataclass, available since python 3.7
    # Always contains the fields "length" and "label".
    # If undefined, length is 0 (zero) and label is an empty string.
    # Different sets of fields may be added during computation.

    def __init__(self, length=0.0, label=""):
        self.length = length
        self.label = label


#############################################################################################
#############################################################################################

class Topostruct():
    """Class that emulates a struct. Keeps topology-related info"""

    # Python note: perhaps replace with dataclass, available since python 3.7
    # Contains the fields "count" and "treestring".
    def __init__(self, count=1, treestring=""):
        self.count = count
        self.treestring = treestring


#############################################################################################
#############################################################################################

class TreeError(Exception):
    pass

#############################################################################################
#############################################################################################

class Tree():
    """Class representing basic phylogenetic tree object."""

    # Implementation note: Tree objects can be constructed from three different kinds of things:
    # (1) Newick tree strings, (2) Bipartition lists, and (3) A list of leaves.
    # It is also possible to create random trees with a given number of leaves.
    # The Tree class therefore has four alternate constructors named:
    # "from_string", "from_biplist", "from_leaves", and "randtree" implemented as classmethods
    # The main constructor "__init__" is therefore mostly empty

    def __init__(self):
        self.parent_dict = None         # Dict node:parent relationships (only built if required)
        self.dist_dict = None
        self.path_dict = None
        self.sorted_intnode_cache = None
        self.bipdictcache = None        # Bipdict cache to avoid multiple calls to bipdict()
        self.topologycache = None       # Topology cache to avoid multiple calls to topology()

    #######################################################################################
    @classmethod
    def from_string(cls, orig_treestring, transdict=None):
        """Constructor 1: Tree object from tree-string in Newick format"""

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
        obj.parent_dict = {}

        # Preprocess parts list to remove any labels or lengths below root node
        # (these are explicitly ignored!)
        # Done by removing all parts between the semicolon and the last right parenthesis
        # NOTE: I am assuming that semicolon is at [-1]. Should I find it first instead?
        while tree_parts_list[-2] != ")":
            del tree_parts_list[-2]

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
                    obj.parent_dict[nodeno] = parent   # Build parent_dict as nodes are added

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
            # If previous part was simultaneously a right parenthesis or a leaf name, then the token
            # must be a branch label
            elif prevpart in [")", "leafname"]:
                child = node_stack[-1]
                parent = node_stack[-2]
                obj.tree[parent][child].label = part

            # Last possibility (I hope): the name is a leafname
            else:
                # If translation dictionary was supplied: change name accordingly
                if transdict:
                    child = transdict[part]
                else:
                    child = sys.intern(part)    # "Intern" strings so only one copy in memory
                                                # (may happen automatically anyway?...)
                    #child = part

                # Check for duplicated leaf names
                if child in obj.leaves:
                    msg = "Leaf name present more than once: {}".format(child)
                    raise TreeError(msg)

                parent=node_stack[-1]                      # Add new leaf node to previous node's
                obj.tree[parent][child] = Branchstruct()   # list of children
                obj.parent_dict[child] = parent            # Build parent_dict as nodes are added

                node_stack.append(child)                   # Push new leaf node onto stack
                obj.leaves.add(child)                      # Also update list of leaves

                part = "leafname"

            prevpart = part

        obj.parent_dict[obj.root] = None

        # If nodes remain on stack after parsing, then something was wrong with tree-string
        if node_stack:
            msg = "Imbalance in tree-string: '%s'" % orig_treestring
            raise TreeError(msg)


        obj.nodes = obj.leaves | obj.intnodes
        return obj

    #######################################################################################

    @classmethod
    def from_biplist(cls, biplist):
        """Constructor 2: Tree object from bipartition list"""

        # Implementation note: @classmethod useful way of implementing alternative constructors

        # Input is a bipartitionlist (actually a dictionary of bipartition:Branchstruct pairs):
        # Names of leaves on one side of a branch are represented as an immutable set of leaves
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree is represented as a dictionary of bipartition:Branchstruct pairs
        obj = cls()

        # Extract set of leaves
        tmpbip = list(biplist.keys())[0]    # This is a set of two leaf name sets (a bipartition)
        part1, part2 = tmpbip               # Get two sets
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
        for (bipart, branchstruct) in biplist.items():
            bip1, bip2 = bipart
            # If this bipartition represents external branch, then update relevant Branchstruct
            if len(bip1) == 1 or len(bip2) == 1:
                if len(bip1) == 1:
                    (leaf, ) = bip1             # A one-member tuple used for value unpacking (pop)
                else:
                    (leaf, ) = bip2

                # Find parent, update branchstruct
                for parent in obj.tree:
                    if leaf in obj.tree[parent]:
                        obj.tree[parent][leaf] = branchstruct
                        break

            # If bipartition represents internal branch: add branch to tree, transfer Branchstruct
            else:
                mrca1 = obj.findMRCA(bip1)
                mrca2 = obj.findMRCA(bip2)

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

    #######################################################################################

    @classmethod
    def from_leaves(cls, leaflist):
        """Constructor 3: star-tree object from list of leaves"""

        treelist = ["("]
        for name in leaflist:
            treelist.append(name)
            treelist.append(",")
        del treelist[-1]
        treelist.append(");")
        return cls.from_string("".join(treelist))

    #######################################################################################

    @classmethod
    def randtree(cls, leaflist=None, ntips=None, randomlen=False, name_prefix="s"):
        """Constructor 4: tree with random topology from list of leaf names OR number of tips"""

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

    #######################################################################################

    def __iter__(self):
        """Returns iterator object for Tree object. Yields subtrees with .basalbranch attribute"""

        # Basal branch struct may contain useful information: e.g., label and length below subtree
        # Single node trees consisting of leaves are also considered subtrees (REMOVE????)
        class Subtree_iterator():
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

        return Subtree_iterator(self)

    #######################################################################################

    def __str__(self):
        """Prints a table of parent-child relationships in the tree including branch lengths and labels"""

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

        # Find widest string in each column (formatting of table could go into utils module...)
        maxwidth = [0]*4
        for i in range(len(table)):
            for j in range(4):
                if len(table[i][j]) > maxwidth[j]:
                    maxwidth[j] = len(table[i][j])

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

    #######################################################################################

    def __eq__(self, other):
        """Implements equality testing for Tree objects"""

        # Two trees are identical if they have the same leaves, the same topology
        # and the same branchlengths (NB: floating point comparison). Branch labels are ignored
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
                if (abs(len1 - len2) / len1) > 0.0000001:       # Floating point comparison
                    return False

        # If we made it this far without returning, then Tree objects must be identical
        return True

    #######################################################################################

    def __hash__(self):
        """Implements hashing for Tree objects, so they can be used as keys in dicts"""

        # Using hash = id of self.
        # NOTE: this does NOT live up to reasonable hash-criteria... Change at some point.
        # NOTE2: Also unsure about effect on performance
        return id(self)

    #######################################################################################

    def build_parent_dict(self):
        """Forces construction of parent_dict enabling faster lookups (avoid function, use dict directly)"""
        self.parent_dict = {}
        for parent in self.intnodes:
            for child in self.tree[parent]:
                self.parent_dict[child] = parent
        self.parent_dict[self.root] = None                  # Add special value "None" as parent of root

    #######################################################################################

    def build_dist_dict(self):
        """Construct dictionary keeping track of all pairwise distances between nodes"""

        # Data structures and algorithm inspired by the Floyd-Warshall algorithm, but modified and faster than O(n^3)
        # since it is on a tree (unique paths)

        # Python note: maybe I could check for existence before recomputing (but then important to clear after changes!)

        # Python note 2: dict.fromkeys does something clever about presizing dict so there is less of a performance
        # hit when it is later added to, hence the slightly odd initialisation below (25% faster than dict comprehension)

        dist = self.dist_dict = dict.fromkeys(self.nodes)
        for key in dist:
            dist[key] = {}
        tree = self.tree
        combinations = itertools.combinations

        # Traverse tree starting from root, breadth-first (sorted_intnodes): ensures below algorithm will be correct
        intnodes = self.sorted_intnodes()
        for parent in intnodes:
            children = tree[parent].keys()
            if dist[parent]:
                prev_contacts = dist[parent].keys()
                for child in children:
                    childlen = tree[parent][child].length
                    for prev_contact in prev_contacts:
                        dist[prev_contact][child] = dist[child][prev_contact] = dist[prev_contact][parent] + childlen
            for (child1, child2) in combinations(children, 2):
                dist[child1][child2] = dist[child2][child1] = tree[parent][child1].length + tree[parent][child2].length
            for child in children:
                dist[parent][child] = dist[child][parent] = tree[parent][child].length

        # Fill in diagonal (zero entries), just in case
        for node in self.nodes:
            dist[node][node] = 0

    #######################################################################################

    def build_path_dict(self):
        """Construct dictionary keeping track of all pairwise paths between nodes"""

        # Data structures and algorithm inspired by the Floyd-Warshall algorithm, but modified and faster than O(n^3)
        # since it is on a tree (unique paths)
        # Python note: dict.fromkeys does something clever about presizing dict so there is less of a performance
        # hit when it is later added to, hence the slightly odd initialisation below (25% faster than dict comprehension)
        path = self.path_dict = dict.fromkeys(self.nodes)
        for key in path:
            path[key] = {}
        tree = self.tree
        combinations = itertools.combinations

        # Traverse tree starting from root, breadth-first (sorted_intnodes): ensures below algorithm will be correct
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

    #######################################################################################

    def sorted_intnodes(self, deepfirst=True):
        """Returns sorted intnode list for breadth-first traversal of tree. Default is to place deep nodes first"""

        # "intnodelist" is a set, meaning iteration occurs in no defined order.
        # This function returns a list sorted such that deep nodes generally go before
        # shallow nodes (deepfirst=False reverses this)

        # If info is in cache, use that. If not: build sorted list and save in cache
        # Python note: maybe use LRU cache instead of keeping track yourself? Then remember to delete when rerooting etc.
        if self.sorted_intnode_cache is not None:
            sorted_nodes = self.sorted_intnode_cache
        else:

            # Add nodes one tree-level at a time. First root, then children of root, then children of those, etc
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

    #####################################################################################

    def leaflist(self):
        """Returns list of leaf names sorted alphabetically"""

        l = list(self.leaves)
        l.sort()
        return l

    #####################################################################################

    def children(self, parent):
        """Returns set containing parent's immediate descendants"""

        # Python note: does not seem to benefit from lru_caching, and leads to multiple problems (I tried)

        try:
            return set(self.tree[parent].keys())
        except KeyError as e:
            msg = "Node %s is not an internal node" % parent
            raise TreeError(msg) from e

    #######################################################################################

    @functools.lru_cache(maxsize=None)
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

    #######################################################################################

    def parent(self, node):
        """Returns parent of node"""

        # First time function is run it will build dictionary of node:parent connections
        # Assumption is that if this function is required once, then it will be used again
        # and building the dict will pay off. Building the dict should not be the default
        # behaviour when new tree is constructed

        # Construct parent_dict if it does not exist already
        if self.parent_dict is  None:
            self.build_parent_dict()

        # Return parent
        try:
            return self.parent_dict[node]
        except KeyError as e:
            raise TreeError("Node {} does not exist".format(node)) from e

    #######################################################################################

    def nearleafs(self, leaf1, maxdist):
        """Returns set of leaves that are less than maxdist away from leaf along tree (patristic distance)"""

        otherleaves = self.leaves - {leaf1}
        neighbors = set()
        for leaf2 in otherleaves:
            if self.nodedist(leaf1, leaf2) < maxdist:
                neighbors.add(leaf2)

        return neighbors

    #######################################################################################

    def nearestNleaves(self, leaf1, N):
        """Returns set of N leaves closest to leaf along tree (patristic distance)"""

        # Python note: numpy.argsort may be faster, but difficult to include ties (n)

        leaflist = self.leaves.copy() - {leaf1}                 # Set of all leaves except leaf1
        leaflist = list(leaflist)
        distlist = self.nodedistlist(leaf1, leaflist)
        dist_leaf_list = list(zip(distlist, leaflist))          # List of (dist, leaf2) tuples
        dist_leaf_list.sort()                                   # Sort on distance (first item in tuple)
        maxdist = dist_leaf_list[N-1][0]                        # Maximum distance to include
        nearest_leaves = set()
        for dist, leaf in dist_leaf_list:                       # Add all leaves with <= maxdist to set (may be more than N)
            if dist <= maxdist:
                nearest_leaves.add(leaf)
            else:
                break

        return nearest_leaves

    #######################################################################################

    def findMRCA(self, leafset):
        """Finds Most Recent Common Ancestor for the provided set of leaves"""

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

        # Build set of all nodes that are ancestral to leafset
        ancestors = set()
        for node in self.intnodes:
            if leafset <= self.remote_children(node):
                ancestors.add(node)

        # Find ancestor with fewest offspring, this is MRCA
        minoffspring = None
        for node in ancestors:
            no_offspring = len(self.remote_children(node))
            if (minoffspring is None) or (no_offspring < minoffspring):
                mrca = node
                minoffspring = no_offspring

        return mrca

    #######################################################################################

    def find_central_leaf(self, leaflist):
        """Finds central leaf for the provided list of leaves.
        Defined as having approximately equal distance to the two farthest leaves in leaflist"""

        nleaves = len(leaflist)
        #  Cluster has one member: return it (yeah, well...)
        if nleaves == 1:
            return leaflist[0]

        #  Cluster has two members: Pick the leaf farthest from root (would the opposite approach be more reasonable?)
        if nleaves == 2:
            if self.nodedist(leaflist[0]) > self.nodedist(leaflist[1]):
                return leaflist[0]
            else:
                return leaflist[1]

        #  Cluster has more than two members: Find two most distant leafs in leaflist (leaves "spreading out" subtree)
        # and then find the leaf that has approximately the same distance to these two (interpreted as being in a sense halfway between them...)
        basenode = self.findMRCA(set(leaflist))
        sub = self.subtree(basenode)
        (d, L1, L2) = sub.diameter(return_leaves=True)

        # Find leaf having the smallest difference between its distance from L1 and L2 (leaf "halfway" between L1 and L2)
        smallest_diff = d       # Pick value certain to be larger than all dist differences
        for leaf in leaflist:
            diff = abs(sub.nodedist(leaf, L1) - sub.nodedist(leaf, L2))
            if diff < smallest_diff:
                central_leaf = leaf
                smallest_diff = diff
        return central_leaf

    #######################################################################################

    def find_common_leaf(self, leaflist):
        """Finds common leaf for the provided list of leaves.
        Defined as having the smallest average distance to remaining leaves (= many close neighbors)."""

        nleaves = len(leaflist)
        #  Cluster has one member: return it (yeah, well...)
        if nleaves == 1:
            return leaflist[0]

        #  Cluster has two members: Pick the leaf farthest from root (would the opposite approach be more reasonable?)
        if nleaves == 2:
            if self.nodedist(leaflist[0]) > self.nodedist(leaflist[1]):
                return leaflist[0]
            else:
                return leaflist[1]

        # Cluster has more than two members: Find leaf having the smallest average distance to remaining leaves
        # This will identify leaves that are part of denser sub-clusters (instead of single outliers)
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

    #######################################################################################

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
                bipart2 = self.remote_children(node2)    # Remotechildren returns set(node) if node is a leaf
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

        msg = "This is not a monophyletic group:\n%s" % leafstring
        raise TreeError(msg)

    #######################################################################################

    @functools.lru_cache(maxsize=None)
    def nodedist(self,node1,node2=None):
        """Returns distance between node1 and node2 along tree (patristic distance)"""

        # Python note: make recursive to gain advantage of keeping intermediate nodedists in cache?

        # Node2 defaults to root if not given
        if node2 is None:
            node2 = self.root

        # Local copies to speed up access
        root = self.root
        self.build_parent_dict()  # Rethink! (in treelib2 always construct?)
        pdict = self.parent_dict
        tree = self.tree

        # Special case when node1 == node2:
        if node1 == node2:
            return 0.0

        # Find path from node1 back to root Keep track of cumulated distances along the way
        child1 = node1
        node1_ancdist = {}
        node1_ancdist[node1] = 0.0
        while child1 != root:
            parent1 = pdict[child1]
            node1_ancdist[parent1] = node1_ancdist[child1] + tree[parent1][child1].length
            child1 = parent1

        # Find path from node2 back to node on node1's path. Keep track of cumulated distances along the way
        child2 = node2
        cumdist2 = 0.0
        while child2 not in node1_ancdist:
            parent2 = pdict[child2]
            cumdist2 += tree[parent2][child2].length
            child2 = parent2

        # Compute combined distance
        nodedist = cumdist2 + node1_ancdist[child2]
        return nodedist

    #######################################################################################

    def nodedistlist(self, node1, nodelist):
        """Returns list of distances between node1 and nodes in nodelist (same order as nodelist)"""

        self.build_dist_dict()
        distlist = []
        for node2 in nodelist:
            distlist.append(self.dist_dict[node1][node2])
        return distlist

    #######################################################################################

    def nodedepth(self,node):
        """Returns depth of node: distance from rightmost leaf-level to node (i.e., depth of rightmost leaf = 0.0) """

        # The depth of node N can be found from the depth of the root as: depth_N = depth_root - dist(root, N)

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

    #######################################################################################

    def nodepath_fromdict(self, node1, node2):
        """Returns path between node1 and node2 along tree, from preconstructed path_dict"""

        try:
            path = [node1]
            n = node1
            while n != node2:
                n = self.path_dict[n][node2]
                path.append(n)
            return path
        except AttributeError as e:
            raise TreeError("The path dictionary has not been constructed. Cannot use nodepath_fromdict") from e

    #######################################################################################

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
        node2path = node2path[:-1]                      # Remove intersection from node2path
        intersect_index = node1path.index(intersect)
        node1path = node1path[:intersect_index+1]       # Remove all that is upstream of intersect in nodepath1

        # Merge paths to form full path
        while node2path:
            lastnode = node2path.pop()
            node1path.append(lastnode)

        # node1path now contains entire path between nodes and is returned
        return node1path

    #######################################################################################

    def length(self):
        """Returns tree length (sum of all branch lengths)"""

        # Python note: Maybe add lru cache, but then remember to clear this cache if leaves are added or removed

        treelength = 0.0
        for node in self.intnodes:
            for child in self.children(node):
                treelength += self.tree[node][child].length

        return treelength

    #######################################################################################

    def height(self):
        """Returns height of tree: Largest root-to-tip distance"""
        nodeset = self.leaves
        node1 = self.root
        maxdist = self.find_most_distant(node1, nodeset)[1]

        return maxdist

    #######################################################################################

    def diameter(self, return_leaves=False):
        """Return diameter: longest leaf-leaf distance along tree. If return_leaves is True: Return tuple with (maxdist, Leaf1, Leaf2)"""
        # Find the two leaves having the largest pairwise distance. Neat, 2-step algorithm:
        # (1) Pick random leaf, L1, find longest path to other leaf, L2
        # (2) Starting at L2, find longest path, this is longest path in tree! (It's true...)

        # Local copy for faster lookup (necessary?)
        nodedist = self.nodedist

        # Step 1: pick random leaf (L1) and find longest path to other leaf (L2)
        L1 = next(iter(self.leaves))       # Get random leaf from set without popping it (this is ugly in python...)
        maxdist = 0.0
        for node2 in self.leaves:
            dist = nodedist(L1, node2)
            if dist > maxdist:
                maxdist = dist
                L2 = node2

        # Step 2: Find longest path starting at L2, this is longest path in tree
        maxdist = 0.0
        for node2 in self.leaves:
            dist = nodedist(L2, node2)
            if dist > maxdist:
                maxdist = dist
                L3 = node2

        # Return requested result
        # Python note: Bad idea to have varying return values. Decide on output and deal with it at other end
        if return_leaves:
            return (maxdist, L2, L3)
        else:
            return maxdist

    #######################################################################################

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

    #######################################################################################

    def average_pairdist(self, leaflist, return_median=False):
        """Return average pairwise distance between leaves in leaflist, measured along tree. Median can be requested"""
        # Python note: better to return list of pair distances, which can then be averaged etc.

        distlist = []
        for leaf1, leaf2 in itertools.combinations(leaflist, 2):
            distlist.append(self.nodedist(leaf1, leaf2))
        if return_median:
            return statistics.median(distlist)
        else:
            return statistics.mean(distlist)

    #######################################################################################

    def average_ancdist(self, leaflist, return_median=False):
        """Return average distance from leaves to their MRCA, measured along tree. Median can be requested"""
        # Python note: better to return list of pair distances, which can then be averaged etc.

        ancnode = self.findMRCA(set(leaflist))
        distlist = []
        for leaf in leaflist:
            distlist.append(self.nodedist(leaf, ancnode))
        if return_median:
            return statistics.median(distlist)
        else:
            return statistics.mean(distlist)

    #######################################################################################

    def shuffle_leaf_names(self):
        """Shuffles the names of all leaves"""

        # NOTE: this is mostly useful for statistical (resampling-style) analysis of whether
        # patterns of clustered leaves are significant

        # Construct list of old names
        oldnames = list(self.leaves)

        # Construct list of new names by shuffling old names
        newnames = list(self.leaves)
        random.shuffle(newnames)

        # Construct intermediate list of names, used to avoid duplicate names during name-changing process
        tmpnames = ["t" + str(i) for i in range(len(oldnames))]

        # Change old names to temporary names
        for (oldname, newname) in zip(oldnames, tmpnames):
            self.rename_leaf(oldname, newname)

        # Change temporary names to new names
        for (oldname, newname) in zip(tmpnames, newnames):
            self.rename_leaf(oldname, newname)

    #######################################################################################

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
        maximal_clades = matching_clades[:]  # Copy of list (can't remove items from list while iterating)
        for clade1 in matching_clades:
            for clade2 in matching_clades:
                if clade1 < clade2:
                    maximal_clades.remove(clade1)
                    break   # No reason to compare clade1 to any more

        return maximal_clades

    #######################################################################################

    def newick(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6):
        """Returns Newick format tree string representation of tree object"""

        # Distances and labels are printed unless user explicitly request no printing
        # NOTE: This could probably be done slightly faster by iteration (instead of recursion)
        # for instance by using a stack, but unlikely that this function will be heavily used...
        def append_children(parentnode):
            """Recursive function that has main responsibility for building Newick tree string"""

            for child in self.children(parentnode):

                dist = self.tree[parentnode][child].length
                label = self.tree[parentnode][child].label

                if child in self.leaves:
                    treelist.append(child)
                    if label and print_leaflabels:
                        treelist.append(" ")
                        treelist.append(label)
                    if printdist:
                        treelist.append(":{num:.{prec}g}".format(num=dist, prec=precision))
                else:
                    treelist.append("(")
                    append_children(child)
                    treelist.append(")")

                    if label and printlabels:
                        treelist.append(label)
                    if printdist:
                        treelist.append(":{num:.{prec}g}".format(num=dist, prec=precision))

                treelist.append(",")

            del treelist[-1]            # Remove last comma when no more siblings

        # EXECUTION STARTS HERE!
        # Newick tree is built left-to-right as list, and finally converted to string
        # Root requires special treatment, rest of tree managed by recursion
        root = self.root
        treelist = ["("]
        append_children(root)
        treelist.append(");")
        treestring = "".join(treelist)
        return treestring

   ########################################################################################

    def nexus(self, printdist=True, printlabels=True, precision=6):
        """Returns nexus format tree as a string"""

        # Construct header
        stringlist = ["#NEXUS\n\nbegin trees;\n   tree nexus_tree = "]

        # Add newick tree string
        stringlist.append(self.newick(printdist, printlabels, precision))

        # Add footer
        stringlist.append("\nend;\n")

        return "".join(stringlist)

   ########################################################################################

    def figtree(self, printdist=True, printlabels=True, precision=6, colorlist=None, color="0000FF"):
        """Returns figtree format tree as a string. Rudimentary - mostly for coloring leaves. Default color=blue"""

        # Implementation note: Rudimentary and mostly for coloring. Should I add figtree section at end with various settings?
        # See saved figtree file for examples

        # Construct header
        stringlist = ["#NEXUS\nbegin taxa\n\tdimensions ntax={};\n\ttaxlabels\n\t".format(len(self.leaflist()))]
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
        stringlist.append(self.newick(printdist, printlabels, precision))

        # Add footer
        stringlist.append("\nend;\n")

        return "".join(stringlist)

    ########################################################################################

    def bipdict(self):
        """Returns tree in the form of a "bipartition dictionary" """

        # Names of leaves on one side of a branch are represented as an immutable set
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree is represented as a dictionary where the keys are bipartitions
        # The values are structs containing the fields "length" (branch length) and
        # "label" (branch label)

        # In the unlikely event this has already been computed: return stored result
        if self.bipdictcache:
            return self.bipdictcache

        bipartition_dict = {}
        leaves = frozenset(self.leaves)

        # For each branch: find bipartition representation, add this and Branchstruct to list.
        # Remote kids of node most distant from root (or node itself) forms one part of bipartition.
        # Other part is then found as diff between all leaves and bipart1
        # Python note: Sorting pays off because remote_children cache is then built in rational order
        sortedintnodes = self.sorted_intnodes(deepfirst=False)

        for node1 in sortedintnodes:
            for node2 in self.children(node1):
                bipart1 = frozenset(self.remote_children(node2))
                bipart2 = leaves - bipart1
                bipartition = frozenset([bipart1, bipart2])

                # Bipartitions are kept in global dictionary to avoid redundant saving (they can be huge).
                # This is especially useful in connection with toposummaries (where potentially tens of
                # thousands of trees all contain tens of bipartitions, but where a large fraction of the
                # bipartitions are used in most trees).
                if bipartition not in Globals.biparts:
                    Globals.biparts[bipartition] = bipartition
                global_bipart = Globals.biparts[bipartition]

                bipartition_dict[global_bipart] = Branchstruct(self.tree[node1][node2].length, self.tree[node1][node2].label)

        # If root is attached to exactly two nodes, then two branches correspond to the same
        # bipartition. Clean up by collapsing two branches (add lengths, compare labels)
        root = self.root
        rootkids = self.children(root)

        if len(rootkids) == 2:

            # First: find out what bipartition root is involved in.
            kid1, kid2 = rootkids
            bipart1 = frozenset(self.remote_children(kid1))
            bipart2 = leaves - bipart1
            bipartition = Globals.biparts[frozenset([bipart1, bipart2])]

            # Overwrite previous info
            bipartition_dict[bipartition] = Branchstruct(self.tree[root][kid1].length, self.tree[root][kid1].label)

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

        self.bipdictcache = bipartition_dict           # Avoid computing several times
        return self.bipdictcache

    ########################################################################################

    def topology(self):
        """Returns set of sets of sets representation of topology ("naked bipartitiondictionary")"""

        # Names of leaves on one side of a branch are represented as an immutable set.
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree topology is represented as a set of bipartitions
        # This is essentially a naked version of a bipdict

        # Is precomputed topology or bipdict present? (Avoid unnecessary computation)
        if self.topologycache:
            return self.topologycache
        elif self.bipdictcache:
            bipdict = self.bipdictcache
        else:
            bipdict = self.bipdict()

        topology = frozenset(bipdict.keys())
        self.topologycache = topology

        return topology

    #######################################################################################

    def is_compatible_with(self, bipart):
        """Checks whether a given bipartition is compatible with the tree"""

        # A bipartition is compatible with a tree if, for all bipartitions in the tree,
        # one half of the query bipart is a subset of one half of the tree bipart.
        query1, query2 = bipart
        topo = self.topology()

        for bipartition in topo:
            tree1, tree2 = bipartition

            # If no query set is a subset of any tree set, then bipart is not compatible with tree
            if not (query1<=tree1 or query1<=tree2 or query2<=tree1 or query2<=tree2):
                return False

        # If we fell of the end of the for loop without returning, then bipart must be compatible
        return True

    #######################################################################################

    def is_resolved(self):
        """Checks whether tree is fully resolved (no polytomies)"""

        n_rootkids = len(self.children(self.root))
        n_leaves = len(self.leaves)
        n_nodes = n_leaves + len(self.intnodes)

        # If root is a bifurcation, and tree resolved, then n_nodes = 2 * n_leaves - 1
        if (n_rootkids == 2) and (n_nodes < (2 * n_leaves - 1)):
            return False

        # If root is a trifurcation ("basal polytomy"), and tree resolved, then n_nodes = 2 * n_leaves - 2
        if (n_rootkids == 3) and (n_nodes < (2 * n_leaves - 2)):
            return False

        # If root has more than 3 kids, then tree is not resolved
        if n_rootkids > 3:
            return False

        # If we made it this far, then tree is resolved
        return True

    #######################################################################################

    def resolve(self):
        """Randomly resolves multifurcating tree by by adding zero-length internal branches."""

        # Find nodes with > 2 children, add to list of nodes needing to be resolved
        unresolved_nodes = []
        # for node in self.sorted_intnodes(deepfirst=False):
        for node in self.intnodes:
            numkids = len(self.tree[node])               # Note: not safe to use .children method while changing tree (cache will break)
            if numkids > 2:
                unresolved_nodes.append(node)

        # Keep adding extra internal nodes until there are no unresolved nodes left
        while unresolved_nodes:
            intnode1 = unresolved_nodes.pop()
            kids = self.tree[intnode1].keys()
            kids = list(kids)

            # Divide children into two random subsets (note: subsets can contain mix of leaves and intnodes)
            subset1_size = random.randint(1, len(kids) - 1)     # Each subset must contain at least 1, at most Nkids - 1
            subset1 = random.sample(kids, subset1_size)
            subset2 = set(kids) - set(subset1)                       # Use set arithmetic to find remaining children
            subset2_size = len(subset2)

            # For each subset:
            #   if more than 1 member: insert branch and intnode2 between intnode1 and subset
            #   if more than 2 members: also add intnode2 to list of unresolved nodes
            # After this, intnode1 is resolved (two branches emanating from it)
            if subset1_size > 1:
                intnode2 = self.insert_node(intnode1, subset1)
                if subset1_size > 2:
                    unresolved_nodes.append(intnode2)

            if subset2_size > 1:
                intnode2 = self.insert_node(intnode1, subset2)
                if subset2_size > 2:
                    unresolved_nodes.append(intnode2)

    #######################################################################################

    def setlength(self, node1, node2, length):
        """Sets length of branch connecting node1 and node2"""

        if node1 == self.parent(node2):
            parent = node1
            child = node2
        elif node2 == self.parent(node1):
            parent = node2
            child = node1
        else:
            msg = "There is no branch connecting node %s and %s" % (node1, node2)
            raise TreeError(msg)

        self.tree[parent][child].length = length

    #######################################################################################

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

    #######################################################################################

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

    #######################################################################################

    def subtree(self, basenode, return_basalbranch=False):
        """Returns subtree rooted at basenode as Tree object. Note: rooting matters! Note 2: basenode may be leaf!"""

        if return_basalbranch:
            if basenode == self.root:
                msg = "Can not return branch below root node"
                raise TreeError(msg)
            parent = self.parent(basenode)
            basalbranch = self.tree[parent][basenode]

        # Special case: basenode is a leaf => subtree is minimal tree with two nodes (root and leaf)
        if basenode in self.leaflist():
            other = Tree.from_string("({});".format(basenode))

        # If basenode is internal: subtree has more than one leaf
        else:
            # Create empty Tree object. Transfer relevant subset of self's data structure to other
            other = Tree()
            other.parent_dict = None
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

    #######################################################################################

    def graft(self, other, node1, node2=None, blen1=0, blen2=0, graftlabel=None):
        """Graft other tree to self. Tree2 (other) intnodes renamed if names clash with those in tree1.
                node1: node in tree1 (self) below which tree2 (other) will be grafted. Must be specified, cannot be root1
                node2: node in tree2 (other) below which tree2 will be attached (defaul is root of tree2)
                blen1: length of branch added to tree1 below graftpoint (lower of two newly created branches)
                blen2: length of branch above graft point and below tree2 (upper of two newly created branches)
                graftlabel: prepend value of "label" to leaf names on t2 (e.g: "graft_s1")"""

        # Add new intnode on branch in tree1 where tree2 is to be grafted
        parent1 = self.parent(node1)
        graftpoint = self.insert_node(parent1, [node1], branchlength=blen1)         # Note: insert_node() requires second argument to be iterable

        # If node2 is not given: set to root node of tree2
        if node2 is None:
            node2 = other.root

        # If node2 is not root of tree2 then re-root on branch below node2. After this, grafting can happen below root2
        elif node2 != other.root:
            other.deroot()
            node2parent = other.parent(node2)
            other.reroot(node2parent, node2)

        # Rename other's internal nodes if names clash with self's
        renameset = self.intnodes & other.intnodes
        if renameset:
            newnum = max(self.intnodes | other.intnodes)        # Start numbering above all numbers used in self and other
            for oldnum in renameset:
                newnum += 1
                other.rename_intnode(oldnum, newnum)

        # prepend label to leaf names on grafted subtree if requested
        if graftlabel is not None:
            for oldname in other.leaflist():
                newname = "{}{}".format(graftlabel, oldname)
                other.rename_leaf(oldname, newname)

        # Update main data structure (self.tree dictionary) by merging with corresponding dictionary from other
        self.tree.update(other.tree)                                        # Merge dictionary from other.tree into self.tree
        self.tree[graftpoint][other.root] = Branchstruct(length=blen2)      # Link subtree to graftpoint in self.tree

        # Update look-up lists and caches
        self.nodes.update( other.nodes )
        self.intnodes.update( other.intnodes )
        self.leaves.update( other.leaves )

        # Reset those caches and lists that are too hard to salvage (or... could I?)
        self.sorted_intnode_cache = None
        self.build_parent_dict()

        # Erase caches
        self.bipdictcache = None
        self.topologycache = None
        self.sorted_intnode_cache = None

        # Clear lru_caches (which cannot be edited manually)
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()


    #######################################################################################

    def cluster_n(self, nclust):
        """Divides tree into 'nclust' clusters based on distance from root.

           Returns tuple containing: list with sets of leafnames (one set per cluster)
                                     list of basenodes of clusters
                                     set containing all leaves that are put in clusters"""

        # Finds clusters by conceptually cutting across branches at height where there are nclust groups downstream of cutpoint
        # This essentially finds equally spaced clusters. Leafs that are upstream of cutpoint (closer to root) also form clusters of their own

        # Parse childlist in .tree, and for each "from" and "to" node, computes distance from root (height)
        # Saves this information as two new atributes on Branchstruct: "parent_height" and "kid_height"

        # Also create sorted list of internal node heights, and number of branches emanating from node as list of (pheight, nkids) tuples
        # This information can be used to infer number of clusters when cutting above a given node in the following way:
        # If two branches emanate from root, then two clusters will be formed by cutting above it. BUT, further out the tree:
        # For each N branches emanating from an internal node, N - 1 additional clusters will be formed by cutting above it
        # This is because one branch goes INTO the node (it was already part of one cluster)
        # (e.g., if two branches emanate from next node futher out, then ONE additional cluster will be formed by cutting downstream of it)

        # Python note: Lazy implementation where I use .nodedist() function for each node in tree. Could probably be sped up
        #               by using information already in .tree dict (thereby essentially looking at lower branches only once)
        # Implementation note: Should this be run each time a new tree is created? And also each time tree is changed then?

        # PYTHON NOTE: simplify. Maybe always return list of leaf-sets. Make sure that root dists are computed from bottom of tree first

        # Sanity check: There can not be more clusters than leaves
        if nclust > len(self.leaves):
            msg = "Error: Requested {} clusters but tree only has {} leaves".format(nclust, len(self.leaves))
            raise TreeError(msg)

        # Find nodeheights and number of emanating branches for all internal nodes
        nodeheightlist = []
        for parent in self.tree:
            pheight = self.nodedist(parent)
            nkids = len(self.tree[parent])
            nodeheightlist.append((pheight, nkids))
            for kid in self.tree[parent]:
                kheight = self.nodedist(kid)
                self.tree[parent][kid].parent_height = pheight
                self.tree[parent][kid].kid_height = kheight
        nodeheightlist.sort()
        self.nodeheightlist = nodeheightlist

        # Find height cutoff corresponding to nclust groups: walk along list of internal nodes until nclust groups reached
        cutno = nodeheightlist[0][1]                    # Initialize cutno to be the number of branches emanating from root
        level = 0
        cutoff = 0.0
        while cutno < nclust:                           # Go out one more node until we have at least required number of clusters
            level += 1
            cutoff = nodeheightlist[level][0]           # Height of node above which we want to cut branches
            cutno += (nodeheightlist[level][1] - 1)     # For each n additional branches cut, there will be n - 1 additional clusters

        # We now know cutoff, i.e., height above which we can cut all branches in order to obtain >= nclust groups
        # PYTHON NOTE: not quite: low leaves (that lie below cutpoint) wont belong to any cluster, so number may be too low. Think!
        # For each branch in tree: Check if cut at cutoff will cross this branch. If so: add node above it to list of base nodes
        clusterlist = []                # list of sets of leaves (each set is one cluster)
        cluster_basenodes = []          # List of basenodes of clusters
        cluster_leaves = set()          # set containing all leaves that are put in clusters
        for parent in self.tree:
            for kid in self.tree[parent]:
                if self.tree[parent][kid].parent_height <= cutoff < self.tree[parent][kid].kid_height:
                    cluster = self.remote_children(kid)
                    clusterlist.append(cluster)
                    cluster_basenodes.append(kid)       # NOTE: some basenodes may be leaves
                    cluster_leaves.update(cluster)
        unclassified = self.leaves - cluster_leaves      # Set containing any unclassified leaves (leaves that are below cutpoint)

        return (clusterlist, cluster_basenodes, unclassified)

    #######################################################################################

    def cluster_cut(self, cutoff):
        """Divides tree into clusters by conceptually cutting across tree at "cutoff" distance from root.
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
        unclassified = self.leaves - cluster_leaves      # Set containing any unclassified leaves (leaves that are below cutpoint)

        return (clusterlist, cluster_basenodes, unclassified)

    #######################################################################################

    def insert_node(self,parent,childnodes,branchlength=0, lab=""):
        """Inserts an extra node between parent and children listed in childnodes list.

        Length of inserted branch is 'branchlength' and defaults to zero. The node number
        of the new node is returned"""

        # Local copies for faster access
        tree = self.tree
        parent_dict = self.parent_dict

        if parent not in self.nodes:
            msg = "Node %d does not exist" % parent
            raise TreeError(msg)

        # Find next un-used node number
        newnode = max(self.intnodes) + 1

        # Add entry for new node in tree
        tree[newnode] = {}

        # Add new internal node as child of "parent"
        tree[parent][newnode] = Branchstruct(length = branchlength, label=lab)

        # Move childnodes from previous parent to new node
        for child in childnodes:
            tree[newnode][child] = tree[parent][child]
            del tree[parent][child]

        # Update self.intnodes and self.nodes to include new node.
        self.intnodes.add(newnode)
        self.nodes.add(newnode)

        # Update self.parent_dict if it exists
        if parent_dict is not None:
            for child in childnodes:
                parent_dict[child] = newnode
            parent_dict[newnode] = parent

        # Erase caches
        self.bipdictcache = None
        self.topologycache = None
        self.sorted_intnode_cache = None

        # Clear lru_caches (which cannot be edited manually)
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()

        return newnode

    #######################################################################################

    def add_branch(self, bipart, blen=0.0, label=""):
        """Adds branch represented by bipartition to unresolved tree."""

        # NOTE: huge overlap with TreeFromBiplist - should use this function there as well!

        # Sanity check: is bipartition compatible with tree?
        if not self.is_compatible_with(bipart):
            raise TreeError("Bipartition is not compatible with tree: %s" % bipart)

        part1, part2 = bipart
        mrca1 = self.findMRCA(part1)
        mrca2 = self.findMRCA(part2)

        # Determine where to insert new node
        # In the special case of a star tree: add two new nodes, moving each half of bipartition
        #       away from root (branch length will be divided equally between two branches)
        if len(self.intnodes) == 1:     # If startree: add two new internal nodes on either side of root
            self.insert_node(self.root, part1, branchlength=blen/2, lab=label)
            self.insert_node(self.root, part2, branchlength=blen/2, lab=label)

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

            # Add branch at determined position: Note - this takes care of updating parent_dict
            self.insert_node(insertpoint, movelist, branchlength=blen, lab=label)

        # Erase caches
        self.bipdictcache = None
        self.topologycache = None

        # Clear lru_caches (which cannot be edited manually)
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()

    #######################################################################################

    def remove_branch(self, node1, node2):
        """Removes branch connecting node1 and node2 (thereby creating polytomy)"""

        # Python note: Length of collapsed branch is lost, so treelength and patristic distances change.
        # Should I distribute lost branch length to retain treelength?
        # Eg, retain treelength and patristic distances within subtree (+ blen/n to each child)
        # Alternatively retain patristic distances to outside (+ blen to each child), but then treelength increases

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
            # Update parent_dict if it exists
            if self.parent_dict is not None:
                self.parent_dict[grandchild] = parent

        # Delete "child" node and link from parent to child. Update intnodes and nodes
        del self.tree[child]
        del self.tree[parent][child]
        self.intnodes.remove(child)
        self.nodes = self.leaves | self.intnodes

        # Update parent_dict if it exists
        if self.parent_dict is not None:
            del self.parent_dict[child]

        # Update self.sorted_intnode_cache if it exists
        if self.sorted_intnode_cache is not None:
            self.sorted_intnode_cache.remove(child)

        # Erase caches
        self.bipdictcache = None
        self.topologycache = None

        # Clear lru_caches (which cannot be edited manually)
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()

    ########################################################################################

    def remove_leaves(self, leaflist):
        for leaf in leaflist:
            self.remove_leaf(leaf)

    ########################################################################################

    def remove_leaf(self, leaf):
        """Removes named leaf from tree, cleans up so remaining tree structure is sane"""

        parent = self.parent(leaf)
        childset = self.children(parent)
        root = self.root
        #print("leaf: {}.   root: {}.   parent: {}.  childset: {}".format(leaf, root, parent, childset)) #DEBUG
        #print("self.intnodes: {}".format(self.intnodes)) #DEBUG
        # If leaf is part of bifurcation AND is directly attached to root, then
        # the "other child" of the root must become the new root
        if (len(childset) == 2) and (leaf in self.children(root)):
            [child2] = childset - {leaf}                    # Remaining item is other child
            del self.tree[root]                             # Remove entry for old root
            self.intnodes.remove(root)
            self.nodes.remove(root)
            self.root = child2                              # child2 is new root
            del self.parent_dict[child2]                    # clean up parent_dict

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

        # Erase caches (possibly something could be salvaged).
        self.bipdictcache = None
        self.topologycache = None

        # Erase self.sorted_intnode_cache if it exists (could I salvage something?)
        self.sorted_intnode_cache = None

        # Clear lru_caches (which cannot be edited manually)
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()

    #######################################################################################

    def collapse_clade(self, leaflist, newname="clade"):
        """Replaces clade (leaves in leaflist) with single leaf. Branch length is set to average dist from basenode parent to leaves"""

        if len(leaflist) == 1:
            # Special case where there is only one leaf in leaflist: Do not collapse anything, but change name to newname (?)
            oldname = leaflist.pop()
            self.rename_leaf(oldname, newname)
        else:
            # Find average distance from parent of basenode to leaves (median - use mean instead?)
            mrca = self.findMRCA(leaflist)
            mrca_parent = self.parent(mrca)
            avdist = self.average_ancdist(leaflist, return_median=True) + self.nodedist(mrca, mrca_parent)

            # Remove all but one of the leaves in leaflist (hackish way of keeping leaf node for subsequent renaming...)
            leaflist = list(leaflist)
            subleaflist = leaflist[1:]
            self.remove_leaves(subleaflist)

            # Rename remaining leaf to newname and set branch length to average dist
            self.rename_leaf(leaflist[0], newname)
            self.setlength(mrca_parent, newname, avdist)

    #######################################################################################

    def nameprune(self, sep="_", keep_pattern=None):
        """Prune leaves based on name redundancy: Find subtrees where all leaves have same start of name (up to first "_")"""

        # A bit of a hack for very specific project... Find way to generalize (e.g. ask for pattern/regular expression to match)
        # Identify namestarts (up to first occurrence of sep) that occur more than once among leaves
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

    #######################################################################################

    def numberprune(self, nkeep, keeplist=None, keep_common_leaves=False, keep_most_distant=False, return_leaves=False, enforceN = False):
        """Prune tree so 'nkeep' leaves remain. Leaves are chosen to be approximately evenly spaced over tree.
        "keeplist" can be used to specify leaves that _must_ be retained.
        'keep_common_leaves' requests preferential retainment of leaves with many neighbors
        (default is to keep leaves that are as equally spaced as possible)
        'keep_most_distant' requests that the two most distant leaves in tree (which spread out the diameter) should be kept
        'return_leaves': return selected leaves, but do not actually prune tree
        'enforceN' enforce exactly N leaves in pruned tree (normally leaves in includelist and most distant are additional to N)"""

        keepset = set()
        if keeplist:
            keepset.update(keeplist)

        # If requested: Find and retain two most distant leaves in tree
        if keep_most_distant:
            (_, L1, L2) = self.diameter(return_leaves=True)
            keepset.update((L1, L2))

        # If enforceN has not been requested: Find N clusters in addition to possible members of keeplist or the two most distant leaves
        if not enforceN:
            clusters = self.cluster_n(nclust=nkeep)[0]          # "clusters" is a list containing sets of leafnames

        # If enforceN: Iteratively find N, N-1, ... clusters until total retained number of leaves (including those in keeplist etc) is == N
        else:
            N = nkeep
            foundN = False
            while not foundN:
                clusters = self.cluster_n(nclust=N)[0]
                coveredclusters = []
                for cluster in clusters:
                    intersection = cluster & keepset
                    if intersection:
                        coveredclusters.append(cluster)
                Nretained = len(keepset) - len(coveredclusters) + len(clusters)
                if Nretained == nkeep:
                    foundN = True
                else:
                    N -= 1

            # We have now identified:   (1) The value of N that will result in the requested number of retained leaves
            #                           (2) The clusters in which members of keepset are located. Remove these clusters before proceeding
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

    #######################################################################################

    def prune_maxlen(self, nkeep, return_leaves=False):
        """Prune tree so the remaining nkeep leaves spread out maximal percentage of branch length (max phylogenetic diversity)"""

        possible_branches = set()     # Possible starting basal branches for next path to leaf (node1 is on path, and node2 is not)
        used_branches = set()         # Branches that are on the path
        keep_leaves = set()           # Leaves to keep in tree

        # Midpoint root to make initialisation simpler
        # (costly - should rewrite algorithm to start anywhere. On the other hand I would need diameter anyway...)
        self.rootmid()

        # Place central data structures and functions in local namespace for faster lookup
        nodedist = self.nodedist
        nodepath = self.nodepath
        remote_children = self.remote_children
        children = self.children

        # Initialise by adding two branches emanating from root to the list of possible starting branches
        # Note: this only works due to midpoint rooting (which ensures root will be on the first longest path between two leafs)
        rootkids = self.children( self.root )
        for child in rootkids:
            possible_branches.add( (self.root, child) )

        # Until we have added nkeep leaves to path: find longest newpath from existing path to leaf, add to path
        while len(keep_leaves) < nkeep:

            # Among possible starting branches: find the one having the max possible distance to a remote child
            maxdist = 0.0
            for (parent, child) in possible_branches:
                for leaf in remote_children(child):
                    if nodedist(parent,leaf) > maxdist:
                        maxdist = nodedist(parent,leaf)
                        n1, n2, keepleaf = parent, child, leaf

            # Add the found leaf to list of leaves. Remove the basal branch that was used from possible starting branches
            keep_leaves.add( keepleaf )
            possible_branches = possible_branches - { (n1, n2) }

            # Update possible_branches and used branches based on newly added path
            newpath = nodepath( n1, keepleaf )
            for i in range( len(newpath) - 1 ):
                parent, child1 = newpath[i], newpath[i+1]
                used_branches.add( ( parent, child1 ) )
                otherkids = children(parent) - {child1}
                for child2 in otherkids:
                    if (parent, child2) not in used_branches:
                        possible_branches.add( (parent, child2) )

        # If requested: return selected leaves without pruning tree
        # Otherwise: prune tree so only leaves in keepset are retained, by removing all others one at a time
        if return_leaves:
            return keep_leaves
        else:
            discard_leaves = self.leaves - keep_leaves
            self.remove_leaves(discard_leaves)

    #######################################################################################

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

    #######################################################################################

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
        self.tree[parent][newname] = self.tree[parent][oldname]
        del self.tree[parent][oldname]

        # Update self.leaves and self.nodes
        self.leaves.add(newname)
        self.leaves.remove(oldname)
        self.nodes.add(newname)
        self.nodes.remove(oldname)

        # Update self.parent_dict if it exists:
        if self.parent_dict is not None:
            self.parent_dict[newname] = self.parent_dict[oldname]
            del self.parent_dict[oldname]

        # All bets are off regarding how the renamed node has altered the other caches
        self.bipdictcache = None
        self.topologycache = None

        # Clear lru_caches (which cannot be edited manually)
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()

    #######################################################################################

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
        if parent is not None:
            self.parent_dict[newnum] = self.parent_dict[oldnum]
            del self.parent_dict[oldnum]
        for child in kidlist:
            self.parent_dict[child] = newnum

        # Clear lru_caches (which cannot be edited manually)
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()

        # All bets are off regarding how the renamed node has altered the other caches
        self.bipdictcache = None
        self.topologycache = None

    #######################################################################################

    def treedist(self, other, normalise=True, verbose=False):
        """Compute symmetric tree distance (Robinson Foulds) between self and other tree. Normalised measure returned by default"""

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
        if verbose:
            return symdif, symdif_norm, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2
        elif normalise:
            return symdif_norm
        else:
            return symdif

    #######################################################################################

    def treesim(self, other, verbose=False):
        """Compute normalised symmetric similarity between self and other tree"""

        symdif, symdif_norm, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2 = self.treedist(other, verbose=True)
        symsim_norm = 1.0 - symdif_norm

        if verbose:
            return symdif, symdif_norm, symsim_norm, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2
        else:
            return symsim_norm

    ########################################################################################

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
            self.bipdictcache = None
            self.topologycache = None
            self.remote_children.cache_clear()
            self.nodedist.cache_clear()
            self.sorted_intnode_cache = None
            self.dist_dict = None

            # Rebuild parent_dict
            self.build_parent_dict()

    ########################################################################################

    def reroot(self, node1, node2=None, polytomy=False, node1dist=0.0):
        """Places new root on branch between node1 and node2, node1dist from node1"""

        # If tree is to be rooted at basal polytomy, then node1 (base of outgroup) is the new root
        if polytomy:
            newroot = node1

        # If polytomy not requested: new root must be inserted on branch between node1 and node2:
        # Determine which node is parent of other, determine branchlength and label for
        # branch between nodes 1 and 2, figure out what the new branch lengths will
        # be after splitting the branch, and finally insert new node. Bail out if the nodes are not neighbors
        else:
            if node2 is None:
                msg = "Need to specify node2 to reroot() method when rooting at bifurcation"
                raise TreeError(msg)
            if node1 == self.parent(node2):
                parent = node1
                child = node2
                parent_to_root_dist = node1dist
                root_to_child_dist = self.tree[parent][child].length - node1dist
            elif node2 == self.parent(node1):
                parent = node2
                child = node1
                parent_to_root_dist = self.tree[parent][child].length - node1dist
                root_to_child_dist = node1dist
            else:
                msg = "Node %s and %s are not neighbors in tree" % (node1, node2)
                raise TreeError(msg)

            newroot = self.insert_node(parent, [child], parent_to_root_dist, self.tree[parent][child].label)
            self.tree[newroot][child].length = root_to_child_dist

        # Things that were already downstream of newroot do not need to be moved,
        # but things that were previously upstream need to be moved downstream
        # which is done by reversing the links on the direct path going back from newroot to oldroot
        oldroot = self.root
        path_to_old = self.nodepath(newroot, oldroot)
        if self.parent_dict is None:
            self.build_parent_dict()
        for i in range( len(path_to_old) - 1 ):
            newparent, oldparent = path_to_old[i], path_to_old[i+1]
            self.tree[newparent][oldparent] = self.tree[oldparent][newparent]
            del self.tree[oldparent][newparent]
            self.parent_dict[oldparent] = newparent

        # Clean up: clear caches, which are now unreliable
        self.bipdictcache = None
        self.topologycache = None
        self.remote_children.cache_clear()
        self.nodedist.cache_clear()
        self.sorted_intnode_cache = None

        # Update root info:
        self.root = newroot

    ########################################################################################

    def rootmid(self):
        """Performs midpoint rooting of tree"""

        # Remove previous root if present at bifurcation
        self.deroot()

        # Sanity check: if tree has zero length, then midpoint rooting is not possible
        if self.length() == 0.0:
            raise TreeError("All branch lengths are zero - midpoint rooting not possible")

        # Find the two leaves having the largest pairwise distance.
        (maxdist, L1, L2) = self.diameter(return_leaves = True)
        midway = maxdist/2.0

        # Get path between L1 and L2
        path = self.nodepath(L1, L2)

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

    ########################################################################################

    def rootminvar(self):
        """Performs minimum variance rooting of tree"""

        # Based on results in:
        # Minimum variance rooting of phylogenetic trees and implications for species tree reconstruction
        # Uyen Mai, Erfan Sayyari, Siavash Mirarab

        # Sanity check: if tree has zero length, then minimum variance rooting is not possible
        if self.length() == 0.0:
            raise TreeError("All branch lengths are zero - minimum variance rooting not possible")

        # Find all pairwise distances between nodes
        # Store all intnode to leaf distances in numpy 2D array (matrix): each row has leaf-dists for one intnode
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
        # Python note: perhaps do this just-in-time when needed (could I omit some computations for leaves for instance?)
        distvar = distmat.var(axis=1)             # axis=1: variance of each row in matrix

        # For each branch in tree: find local minimum variance point using equation 7 in Mai paper (mimimum of parabola)
        # Location of point is stored as triplet: parentnode, childnode, distance of point from parent on branch
        # Keep track of overall minimum variance and corresponding location

        # Arbitrarily initialise overall minimum to be distance zero from first intnode on one of its branches
        minvar = distvar[0]
        minparent = nodes[0]
        minchild = self.children(minparent).pop()
        minpardist = 0.0

        # Iterate over branches in tree, compute minimum variance point, and update overall minimum if relevant
        # Names of intermediate variables set to match those in Mai paper
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
                x = -b / (2*a)                  # Equation 4. x = distance from parent on branch
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

        # If minparent is not root: remove old root and reroot using the minparent and minpardist already found:
        if minparent != self.root:
            self.deroot()
            self.reroot(node1=minparent, node2=minchild, node1dist=minpardist)

        # If minparent IS root: remove old root and reroot using the minparent and minpardist already found:
        else:
            # Minpardst == 0: Do nothing. Current root is minimal variance root
            if minpardist == 0.0:
                return
            # Minpardist != 0: Do something, depending on whether root is at bifurcation or not
            else:
                rootkids = self.children(self.root)
                # Bifurcation: old root will be removed. Find new minparent (one of root's other children)
                if len(rootkids) == 2:
                    rootkids.remove(minchild)
                    minparent = rootkids.pop()      # Chose (one of) root's other children as minparent
                    minpardist = minpardist + self.dist_dict[self.root][minparent]
                    self.deroot()
                    self.reroot(node1=minparent, node2=minchild, node1dist=minpardist)
                    return

                # Multifurcation: old root will not be removed. Do not need to find new minparent.
                else:
                    self.reroot(node1=minparent, node2=minchild, node1dist=minpardist)
                    return


    ########################################################################################

    def rootout(self,outgroup, polytomy=False):
        """Roots tree on outgroup"""

        # Remove previous root if present at bifurcation
        self.deroot()

        # Find pair of internal nodes corresponding to ingroup:outgroup bipartition
        if type(outgroup) is str:       # if just a single name is given: protect string from iteration to single letters
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

    ########################################################################################

    def spr(self,subtree_node, regraft_node):
        """Subtree Pruning and Regrafting.

            subtree_node: basenode of subtree that will be pruned.
            regraft_node: node in tree below which subtree will be grafted. Must be specified, cannot be root"""

        # Subtree Prune: Create Tree object corresponding to subtree. Remove subtree from self (one leaf at a time)
        subtree = self.subtree(subtree_node)
        for leaf in self.remote_children(subtree_node):
            self.remove_leaf(leaf)

        # Regraft: Add subtree back onto remaining tree
        self.graft(subtree, regraft_node)

#############################################################################################
#############################################################################################
#############################################################################################

class Tree_set():
    """Class for storing and manipulating a number of trees"""

    def __init__(self):
        self.treelist = []

    ########################################################################################

    def __getitem__(self, index):
        """Implements indexing of treeset. Simple index returns single tree. Slice returns Tree_set object with selected subset of trees"""
        if isinstance(index, slice):
            newtreeset = Tree_set()
            newtreeset.treelist = self.treelist[index]
            return newtreeset
        else:
            return self.treelist[index]

    ########################################################################################

    def __len__(self):
        return len(self.treelist)

    ########################################################################################

    def addtree(self, tree):
        """Adds Tree object to Treeset object"""
        self.treelist.append(tree)

    ########################################################################################

    def addtreeset(self, treeset):
        """Adds all trees in Tree_set object to this Tree_set object"""
        for tree in treeset:
            self.addtree(tree)

    ########################################################################################

    def rootmid(self):
        """Performs midpoint rooting on all trees in Tree_set"""
        for tree in self.treelist:
            tree.rootmid()

    ########################################################################################

    def nexus(self, printdist=True, printlabels=True):
        """Returns nexus format tree as a string"""

        # Construct header
        stringlist = ["#NEXUS\n\nbegin trees;\n"]

        # Add newick tree strings
        for tree in self.treelist:
            stringlist.append("    tree nexus_tree = ")
            stringlist.append(tree.newick(printdist, printlabels))
            stringlist.append("\n")

        # Add footer
        stringlist.append("end;\n")

        return "".join(stringlist)

    ########################################################################################

    def newick(self, printdist=True, printlabels=True):
        """Returns newick format tree as a string"""

        stringlist = []
        for tree in self.treelist:
            stringlist.append(tree.newick(printdist, printlabels))
            stringlist.append("\n")
        return "".join(stringlist)


#############################################################################################
#############################################################################################
#############################################################################################

class SmallTreeSummary():
    """Class summarizing bipartitions and branch lengths (but not topologies) from many trees"""

    def __init__(self, include_zeroterms=False):
        """TreeSummary constructor. Initializes relevant data structures"""

        self.tree_count = 0
        self.tree_weight_sum = 0.0
        self.bipartsummary = {}                 # Main repository for bipartition info
        self.bipart_cache = {}                  # Cache storing most frequent bipartitions
        self.cache_minimum = 0.04               # Minimum frequency for bipart to be included in above cache
                                                # it is assumed that allcompat tree can be built solely with
                                                # bipartitions that are more frequent than this!!!!

        self.bipart_processed = False               # Flag indicating whether bipartsummary has been
                                                    # processed (i.e., freq+var has been computed. and
                                                    # bipart_cache has been constructed)

        self.include_zeroterms = include_zeroterms  # Flag indicating whether to count absent branches

    ########################################################################################

    def add_tree(self, curtree, weight=1.0):
        """Add tree object to treesummary, update all relevant summaries"""

        # Unset flag if it has previously been set (happens when adding more trees after
        # calling bipart_result() already)
        self.bipart_processed = False

        # Main interface to TreeSummary.
        # Takes tree object, updates relevant measures
        # First time entered: build set of leaves for consistency checking:
        if self.tree_count == 0:
            self.leaves = curtree.leaves
        elif curtree.leaves != self.leaves:
            msg = "Tree number %d has different set of leaves than previous trees" % self.tree_count
            raise TreeError(msg)

        self.tree_count += 1
        self.tree_weight_sum += weight       # The weighted equivalent of tree_count
        bipdict = curtree.bipdict()

        # I am interested in being able to compute weighted frequency of a bipartition as well as
        # the weighted mean and weighted variance of the branch length for that bipartition.
        # In order to do this I follow the robust (= no underflow/overflow problems), one-pass approach described
        # in D.H.D. West, "Updating Mean and Variance Estimates: An Improved Method", Communications of the ACM,
        # 22(9), 1979. A number of variables are used, some of which are stored in the summary dictionary.
        # These variables have mostly been named according to that paper. (exceptions are "bip_count" which was "N",
        # "weight" which was "W", "mean" which was "M", and "brlen" which was "X").
        for bipart in bipdict:
            brlen = bipdict[bipart].length

            # If bipartition has been seen before: update existing info
            if bipart in self.bipartsummary:
                Q = brlen - self.bipartsummary[bipart].mean
                TEMP = self.bipartsummary[bipart].SUMW + weight
                R = Q*weight/TEMP
                self.bipartsummary[bipart].mean += R
                self.bipartsummary[bipart].T += R*self.bipartsummary[bipart].SUMW*Q
                self.bipartsummary[bipart].SUMW = TEMP
                self.bipartsummary[bipart].bip_count += 1

            # If bipartition has never been seen before: add it to dict and enter info
            else:
                self.bipartsummary[bipart]=Branchstruct()
                self.bipartsummary[bipart].bip_count = 1
                self.bipartsummary[bipart].SUMW = weight
                self.bipartsummary[bipart].mean = brlen
                self.bipartsummary[bipart].T = 0.0



    ########################################################################################

    def update(self, treesummary):
        """Merge this class with external treesummary"""

        # Sanity check: do two treesummaries refer to same set of leaves?
        if self.leaves != treesummary.leaves:
            msg = "Not all trees have same set of leaves."
            raise TreeError(msg)

        # Unset flag if it has previously been set (relevant when adding treesummary after
        # previously calling bipartResult)
        self.bipart_processed = False

        # Update treecount and weight
        self.tree_count += treesummary.tree_count
        self.tree_weight_sum += treesummary.tree_weight_sum

        # Merge "treesummary.bipartsummary" with "self.bipartsummary"
        bisum = treesummary.bipartsummary
        for bipart in bisum:
            # If bipart already in self.bipartsummary, update fields
            if bipart in self.bipartsummary:

                sumw1 = self.bipartsummary[bipart].SUMW
                sumw2 = bisum[bipart].SUMW
                mean1 = self.bipartsummary[bipart].mean
                mean2 = bisum[bipart].mean
                t1 = self.bipartsummary[bipart].T
                t2 = bisum[bipart].T

                self.bipartsummary[bipart].bip_count += bisum[bipart].bip_count
                self.bipartsummary[bipart].mean = (mean1*sumw1 + mean2*sumw2)/(sumw1+sumw2)
                self.bipartsummary[bipart].SUMW += bisum[bipart].SUMW

                # Note: following expression arrived at empirically!
                # I have not proven this is correct, but it does seem to work...
                self.bipartsummary[bipart].T = t1 + t2 + sumw1*sumw2*(mean2-mean1)*(mean2-mean1)/(sumw1+sumw2)

            # If bipartition has never been seen before, simply transfer data from treesummary.bipartsummary:
            else:
                self.bipartsummary[bipart]=Branchstruct()
                self.bipartsummary[bipart].SUMW = bisum[bipart].SUMW
                self.bipartsummary[bipart].mean = bisum[bipart].mean
                self.bipartsummary[bipart].T = bisum[bipart].T
                self.bipartsummary[bipart].bip_count = bisum[bipart].bip_count

    ########################################################################################

    def bipart_result(self):
        """Return raw summary of all observed bipartitions"""

        # Results are stored in the dictionary built in add_tree. Keys are bipartitions,
        # values are structs (classes) that contain a number of intermediate values.
        # Here I use these values to compute and add the fields "freq" (bipartition frequency),
        # and "sem" (standard error of the mean). "mean" (mean branch length) is already present
        # NOTE: I don't actually check if branch lengths are present. If not then mean=0 and sem=0
        if self.tree_count == 0:
            msg = "No trees added to summary object. Impossible to compute result."
            raise TreeError(msg)

        # Compute (1) branch freq, (2) standard error of the mean of branch length.
        for bipart in self.bipartsummary:
            self.bipartsummary[bipart].freq = self.bipartsummary[bipart].SUMW/self.tree_weight_sum

            # If "include_zeroterms" is set, then branch length is considered to be zero in those trees
            # where the bipartition is absent, and these terms are included in the computation
            # (this is done by using the count and weight-sum from all trees instead of from given bipartition)
            if self.include_zeroterms:
                n = self.tree_count
                sumw = self.tree_weight_sum
                # Correct the already computed mean, to account for missing zero terms
                self.bipartsummary[bipart].mean = self.bipartsummary[bipart].mean*self.bipartsummary[bipart].SUMW/sumw
            else:
                n = self.bipartsummary[bipart].bip_count
                sumw = self.bipartsummary[bipart].SUMW
            T = self.bipartsummary[bipart].T

            # The variance (and standard error of the mean) is only defined if n>1.
            # If n==1 var and sem are arbitrarily set to 999999
            if n>1:
                variance = T*n/((n-1)*sumw)         # Weighted variance (ordinary variance if all w=1)
                self.bipartsummary[bipart].var = variance
                self.bipartsummary[bipart].sem = math.sqrt(variance)/math.sqrt(n)    # Standard error of the mean
            else:
                self.bipartsummary[bipart].var = 999999
                self.bipartsummary[bipart].sem = 999999

            # Save bipartition to cache if its frequency is>cache_minimum
            # This cache is the basis for constructing consensus trees, and it is therefore assumed
            # that an allcompat tree can be fully resolved using only bipartitions in this cache.
            # If cache_minimum is 0.01 then this is probably almost always true, but I don't know!!!!
            # Only relevant fields are saved
            if self.bipartsummary[bipart].freq > self.cache_minimum:
                self.bipart_cache[bipart]=Branchstruct()
                self.bipart_cache[bipart].mean = self.bipartsummary[bipart].mean
                self.bipart_cache[bipart].var = self.bipartsummary[bipart].var
                self.bipart_cache[bipart].sem = self.bipartsummary[bipart].sem
                self.bipart_cache[bipart].freq = self.bipartsummary[bipart].freq

        # Set flag indicating that bipartsummary has been processed
        # (meaning that freq+var has been computed, and that bipart_cache has been constructed)
        self.bipart_processed = True

        return self.bipartsummary

    ########################################################################################

    def bipart_to_string(self, bipartition, position_dict):
        """Takes bipartition (set of two leaf sets) and returns string representation"""

        # Meant to be used only by bipartReport. Make pseudo-private?

        bipart1, bipart2 = bipartition

        # Bipartstring will be built as a list of chars first. Initialize with all "."
        stringwidth = len(self.leaves)
        bipart_list = stringwidth * ["."]

        # Smaller set is represented by "*" characters. Larger set by "." characters (already set)
        if len(bipart1) < len(bipart2):
            smallset = bipart1
        else:
            smallset = bipart2

        for leaf in smallset:
            pos = position_dict[leaf]
            bipart_list[pos] = "*"

        return "".join(bipart_list)        # Concatenate into one string

    ########################################################################################

    def bipart_report(self, includeleaves=True, minfreq=0.05):
        """Return processed, almost directly printable, summary of all observed bipartitions"""

        # Bipart report consists of a tuple containing:
        #       (0) a sorted list of leaves (for interpreting bipartstring)
        #       (1) a sorted list of lists. Each item list is: [bipartstring, freq, mean, var, sem]
        #           entire list is sorted on bipartition frequency
        # The boolean argument "includeleaves" controls whether or not to report external branches

        # Compute raw results if they are not already stored
        if not self.bipart_processed:
            self.bipart_result()

        # If minfreq >= cache_mimimum, then all required info is in the much smaller cache - use that!
        if minfreq >= self.cache_minimum:
            raw_result = self.bipart_cache
        else:
            raw_result = self.bipartsummary

        # Must first figure out which leaves correspond to which positions in bipartstring
        # Note: leaves are ordered alphabetically, meaning first char in bipstring corresponds
        # to first leaf in alphabetic sort
        leaflist = sorted(self.leaves)

        position_dict = {}
        for i in range(len(leaflist)):
            leaf = leaflist[i]
            position = i
            position_dict[leaf] = position

        # Loop over all bipartitions in raw_result, build formatted result list in process
        bipreport = []
        for bipart in raw_result:
            bipstring = self.bipart_to_string(bipart, position_dict)
            bipsize = bipstring.count("*")              # Size of smaller set

            # Only report bipartitions that occur more often than "minfreq":
            if raw_result[bipart].freq>minfreq:

                # Only include external branches if "includeleaves" is set
                if bipsize>1 or includeleaves:
                    freq = raw_result[bipart].freq
                    mean = raw_result[bipart].mean
                    var = raw_result[bipart].var
                    sem = raw_result[bipart].sem
                    bipreport.append([freq, bipstring, mean, var, sem])

        # Sort bipreport according to (1) frequency (higher values first), (2) size of
        # smaller bipartition (external branches before internal branches), and
        # (3) bipartstring (*.. before .*. before ..*)
        # First construct temporary list of (1-freq, bipsize, bipstring, originial list-item)
        # tuples. Sort this list of tuples and then re-extract the original list-item again
        # (Example of Decorate, Sort, Undecorate idiom)
        tmplist = sorted([(1-bip[0], bip[1].count("*"), bip[1], bip) for bip in bipreport])
        bipreport = [tup[-1] for tup in tmplist]        # Last element of tuple is orig list

        # Return tuple of (leaflist, bipreport)
        return (leaflist, bipreport)

    ########################################################################################

    def contree(self, cutoff=0.5, allcompat=False, lab = "freq"):
        """Returns a consensus tree built from selected bipartitions"""

        # Check validity of "lab" argument
        if lab not in ["freq", "sem", "rse"]:
            msg = "contree method called with invalid lab argument: %s.\nMust be 'freq', 'sem', or 'rse'" % lab
            raise TreeError(msg)

        if cutoff < 0.5:
            msg = "Consensus tree cutoff has to be at least 0.5"
            raise TreeError(msg)

        # Construct dictionary of most frequent bipartitions.
        # Use bipart_cache if it already exists, if not then construct it first
        if self.bipart_processed:
            bipdict = copy.deepcopy(self.bipart_cache)
        else:
            self.bipart_result()
            bipdict = copy.deepcopy(self.bipart_cache)

        # Initialize new bipdict for keeping relevant bipartitions
        conbipdict = {}

        # Iterate over all bipartitions. Copy info for biparts where freq>cutoff
        # PYTHON NOTE: I am changing bipdict during this iteration, therefore I have to use
        # "for bipart in bipdict.keys()" (which uses a full copy of the keys as a list) and not
        # the slightly faster "for bipart in bipdict:" (which uses the builtin iterkeys method).
        for bipart in list(bipdict.keys()):

            if bipdict[bipart].freq > cutoff:
                # Option "lab" indicates what to use as branch labels. The possible accepted values are:
                #   "freq": bipartition frequency  [DEFAULT]
                #   "sem": standard error of the mean for the branchlength
                #   "rse": relative standard error = sem/mean
                if lab == "freq":
                    conbipdict[bipart] = Branchstruct(length=bipdict[bipart].mean,
                                            label= "%.2f" % bipdict[bipart].freq)
                elif lab == "sem":
                    conbipdict[bipart] = Branchstruct(length=bipdict[bipart].mean,
                                            label= "%.6f" % bipdict[bipart].sem)
                elif lab == "rse":
                    # If branch length is zero: set rse to be "NaN"
                    if bipdict[bipart].mean == 0.0:
                        conbipdict[bipart] = Branchstruct(length=bipdict[bipart].mean,
                                            label= "NaN")
                    # If branch length > 0: set rse = sem/len
                    else:
                        conbipdict[bipart] = Branchstruct(length=bipdict[bipart].mean,
                                            label= "%.2f" % (bipdict[bipart].sem/bipdict[bipart].mean))
                del bipdict[bipart]

        # Build tree from bipartitions  in new bipdict
        contree = Tree.from_biplist(conbipdict)

        # If allcompat has been requested: add remaining, compatible bipartitions to contree
        if allcompat:

            # Construct sorted list of tuples: (frequency, bipart)
            freqbiplist = sorted([(bipdict[b].freq, b) for b in bipdict])
            freqbiplist.reverse()

            # Iterate over sorted list, adding most frequent, compatible bipartitions first
            for _, bipart in freqbiplist:
                if contree.is_compatible_with(bipart):
                    blen = bipdict[bipart].mean
                    if lab == "freq":
                        label= "%.2f" % bipdict[bipart].freq
                    elif lab == "sem":
                        label= "%.6f" % bipdict[bipart].sem
                    elif lab == "rse":
                        label = "%.2f" % bipdict[bipart].freq
                    contree.add_branch(bipart, blen, label)

        return contree


########################################################################################
########################################################################################
########################################################################################

class BigTreeSummary(SmallTreeSummary):
    """Class summarizing bipartitions, branch lengths, and topologies from many trees"""

    # Does everything SmallTreeSummary does and also keeps track of topologies
    # (topology list is potentially quite big, which is the reason for not including it in STS)

    def __init__(self, include_zeroterms=False, outgroup=None, rootmid=False):
        """TreeSummary constructor. Initializes relevant data structures"""

        # Most stuff done by superclass constructor
        SmallTreeSummary.__init__(self, include_zeroterms)

        # This is where topology information is kept
        self.toposummary = {}

        # Save outgroup for rooting (or None if outgroup rooting not requested)
        self.outgroup = outgroup

        # Save flag for whether to do midpoint rooting
        # IMPLEMENTATION NOTE: should check for clash with outgroup rooting
        self.rootmid = rootmid

    ########################################################################################

    def add_tree(self,curtree, weight=1.0):
        """Add tree to treesummary, update all summaries"""

        # Superclass method takes care of updating n_trees and all bipart-related info
        SmallTreeSummary.add_tree(self, curtree, weight)

        # If topology has never been seen before, then add it and initialize count
        # If topology HAS been seen before then update count
        topology = curtree.topology()
        if topology in self.toposummary:
            self.toposummary[topology].count += weight
        else:
            # Attempt to root on specified outgroup, but use midpoint rooting if this fails
            # NOTE: should I warn user?
            if self.outgroup:
                try:
                    curtree.rootout(self.outgroup)
                except TreeError:
                    curtree.rootmid()

            elif self.rootmid:
                curtree.rootmid()
            treestring = curtree.newick(printdist=False, printlabels=False)
            self.toposummary[topology]=Topostruct(weight, treestring)


    ########################################################################################

    def update(self, treesummary):
        """Merge this class with other treesummary"""

        # Superclass method takes care of updating tree_count, tree_weight_sum, and
        # bipartsummary, and also unsets self.bipart_processed flag
        SmallTreeSummary.update(self, treesummary)

        # Merge "treesummary.toposummary" with "self.toposummary"
        topsum = treesummary.toposummary
        for topology in topsum:
            # If topology already in self.toposummary, update count
            if topology in self.toposummary:
                self.toposummary[topology].count += topsum[topology].count
            # If topology has never been seen before, simply transfer entry
            else:
                self.toposummary[topology]=topsum[topology]

    ########################################################################################

    def topo_report(self):
        """Returns list of [freq, treestring] lists"""

        if self.tree_count == 0:
            msg = "No trees added to summary object. Impossible to compute result."
            raise TreeError(msg)

        toporeport = []
        for topology in self.toposummary:
            freq = self.toposummary[topology].count/self.tree_weight_sum
            treestring = self.toposummary[topology].treestring

            # Add freq and treestring to current list
            toporeport.append([freq, treestring])

##            # Delete entry from toposummary to save memory
##            del self.toposummary[topology]

        # Sort report according to frequency (higher values first) and return
        # (In Python, sequence comparison always looks at first element (freq) first)
        toporeport.sort()
        toporeport.reverse()

##        del self.toposummary
        return toporeport


########################################################################################
########################################################################################
########################################################################################


class Treefile():
    """Abstract base-class for representing tree file objects."""

    # Classes for specific formats inherit from this class and add extra stuff as needed.
    # NOTE: i am opening files in "read text" with encoding UTF-8. Will this work across platforms?

    def __init__(self, filename=None, ishandle=False, data=None):

        # Special filename "-" indicates stdin.
        # Read all data from stdin, place it in virtual (RAM) file, and use that instead of real file
        if ishandle:
            self.treefile = StringIO(data)
        elif filename == "-":
            data = sys.stdin.read()
            self.treefile = StringIO(data)
        else:
            self.treefile = open(filename, mode="rt", encoding="UTF-8")

        self.buffer = ""                # Used for keeping leftovers after reading whole line

    ########################################################################################

    def get_treestring(self):
        """Return next tree-string"""

        # We are now at beginning of treestring, read until semi-colon encountered
        # Raise StopIteration when EOF has been reached
        stringlist = [self.buffer]

        if ";" not in self.buffer:           # Only read on if end of treestring not reached
            for line in self.treefile:
                stringlist.append(line)
                if ";" in line: break

        treestring = "".join(stringlist)

        # If we got this far and still haven't found ";" then EOF must have been reached
        # Signal EOF to caller, which can then clean up and stop iteration
        if ";" not in treestring:
            return None

        stringparts = treestring.split(";")
        treestring = "".join([stringparts[0], ";"])
        self.buffer = "".join(stringparts[1:])

        return treestring

    ########################################################################################

    def read_trees(self, discardprop=0.0):
        """Reads trees from file and returns as Tree_set object. Can discard fraction of trees"""

        # Avoid erroneous error message when running pylint ("self" is OK for iteration)
        # pylint: disable=not-an-iterable

        treeset = Tree_set()
        for tree in self:
            treeset.addtree(tree)
        if discardprop > 0:
            ndiscard = int(discardprop * len(treeset))
            treeset = treeset[ndiscard:]

        return treeset

########################################################################################
########################################################################################
########################################################################################


class Newicktreefile(Treefile):
    """Class representing Newick tree file. Iteration returns tree-objects"""

    def __init__(self, filename, ishandle=False, data=None):
        Treefile.__init__(self, filename, ishandle, data)
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

    ########################################################################################

    def __iter__(self):
        return self

    ########################################################################################

    def __next__(self):
        treestring = self.get_treestring()
        if treestring is None:
            self.treefile.close()
            raise StopIteration
        else:
            treestring =  remove_comments(treestring, leftdelim = "[", rightdelim="]")
            return Tree.from_string(treestring)

########################################################################################
########################################################################################
########################################################################################


class Nexustreefile(Treefile):
    """Class representing Nexus tree file. Iteration returns tree object or None"""

    ########################################################################################

    def __init__(self, filename):
        """Read past NEXUS file header, parse translate block if present"""

        Treefile.__init__(self, filename)

        # Can be called with a file-object or any other object that supports
        # iteration by line ("for line in object:") while retaining state information
        # (a second for loop should start where the first for-loop stopped iterating)
        # Optional argument "noreturn" is a boolean that controls whether iteration
        # should return a tree object (default, happens when noreturn==False),
        # or whether None should be returned (happens when noreturn==True)
        # Useful for skipping part of treefile as quickly as possible

        ####################################################################################

        def skip_comment(line):
            """Reads past a NEXUS comment"""

            if line.count("[") > line.count("]"):
                for extraline in self.treefile:
                    line += extraline
                    if line.count("[") == line.count("]"): break

            # Now we should have an equal number of "[" and "]"
            # Remove comments
            return remove_comments(line, leftdelim = "[", rightdelim="]")

        ####################################################################################

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

        ####################################################################################

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
            pattern = re.compile("^.*translate\s*", re.IGNORECASE | re.DOTALL)
            self.buffer = pattern.sub("", self.buffer)

            # Remove end of buffer, only "translate" block core remains after this
            # NOTE: originally had following complicated pattern, which was found to break on some trees.
            # Why did I do that? Are there some cases that I am no longer covering?
            # pattern = re.compile(";\s+u?tree\s+(\*\s)?\s*[\w\-\/\.]+\s*=.*", re.IGNORECASE | re.DOTALL)
            pattern = re.compile(";\s+tree.*", re.IGNORECASE | re.DOTALL)
            transblock = pattern.sub("", self.buffer)

            # Split on commas, generating list of "code-whitespace-origname" pairs
            # then split each of these on whitespace, and build translation dictionary
            translist = transblock.split(",")
            for transpair in translist:
                (code, realname) = transpair.split()
                self.transdict[code] = realname

        # Keep everything after "tree <NAME> =" in self.buffer: may contain tree!
        pattern = re.compile("^.*?=[^(]*", re.DOTALL)        # Non-greedy mathching *?: find first = and then every non-start parenthesis
        self.buffer = pattern.sub("", self.buffer)

    ########################################################################################

    def __iter__(self):
        return self

    ########################################################################################

    def __next__(self, noreturn=False):

        treestring = self.get_treestring()
        if treestring is None:
            self.treefile.close()
            raise StopIteration
        # remove comments in brackets if present
        # remove leading "tree NAME =" (compiled regexp "tree_header_pattern")
        treestring = remove_comments(treestring, leftdelim = "[", rightdelim="]")   #DEBUG: does not deal with figtree comments!
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

########################################################################################
########################################################################################
########################################################################################
##
##          Treebuilding section: classes and methods for building trees
##
########################################################################################
########################################################################################

class Dist_tree():
    """Class for constructing trees using distance based methods"""

    # Can be constructed from either alignment or distance matrix, so class has two
    # alternate constructors
    def __init__(self):
        pass

    #######################################################################################

    @classmethod
    def from_alignment(cls, alignment, dist="pdist"):
        """Constructor 1: constructs distance tree constructor from alignment object"""
        self = cls()
        self.distmat = alignment.distmatrix(dist)
        self.leaves = alignment.getnames()
        return self

    #######################################################################################

    @classmethod
    def from_distmat(cls, distmatrix):
        """Constructor 2: constructs distance tree constructor from distance matrix object"""
        self = cls()
        self.distmat = distmatrix
        self.leaves = self.distmat.getnames()
        return self

    #######################################################################################

    def nj(self):
        """Computes neighbor joining tree, returns Tree object"""

        # Implementation node: Maybe this should be written to resemble UPGMA method below,
        # i.e., place responsibility for merging distances etc in distmatrix object, instead of here?

        # Construct star-tree. This will be resolved node by node during algorithm
        njtree = Tree.from_leaves(self.leaves)
        rootnode = njtree.root

        # Local copy of leaflist for keeping track of (as yet) unclustered nodes
        remaining_nodes = self.leaves[:]

        # Compute "udists": summed dist to all other nodes
        nodes = list(self.distmat.getnames())
        nnodes = len(nodes)
        udist = dict.fromkeys(nodes, 0.0)  # Initialize dict: keys are leaves, vals are 0.0
        for i in range(nnodes):
            n1 = nodes[i]
            for j in range(i+1, nnodes):
                n2 = nodes[j]
                dist = self.distmat.getdist(n1, n2)
                udist[n1] += dist
                udist[n2] += dist

        # Avoid dot-notation to speed up execution
        distmat = self.distmat
        d = self.distmat.getdist
        u = udist

        # Main loop: continue merging nearest neighbors until only two nodes left
        while len(remaining_nodes) > 2:

            nnodes = len(remaining_nodes)
            n_minus_2 = nnodes - 2                   # Move computation out of loop to save time
            smallest = None

            # Find nearest neighbors according to nj-dist: njdist = (n-2) * d(n1, n2) - u(n1) - u(n2)
            for i in range(nnodes):
                n1 = remaining_nodes[i]
                u1 = u[n1]                          # Move lookup out of loop to save time
                for j in range(i+1, nnodes):
                    n2 = remaining_nodes[j]
                    njdist = n_minus_2 * d(n1, n2) - u1 - u[n2]

                    # If "smallest" is not defined: set to current values
                    # If "smallest" is defined and njdist < smallest: update values
                    if (not smallest) or (njdist < smallest):
                        smallest = njdist
                        nb1, nb2 = n1, n2

            # Connect two nearest nodes

            # (1) Update tree, compute length of branches from new node to merged nodes
            newnode = njtree.insert_node(rootnode, [nb1, nb2])
            dist1 = 0.5 * d(nb1, nb2) + 0.5 * (u[nb1] - u[nb2]) / (nnodes - 2)
            dist2 = 0.5 * d(nb1, nb2) + 0.5 * (u[nb2] - u[nb1]) / (nnodes - 2)
            njtree.setlength(newnode, nb1, dist1)
            njtree.setlength(newnode, nb2, dist2)

            # (2) Update distance matrix and list of remaining nodes
            # Python note: I am changing distance matrix here - should I use copy fo I can re-use?
            remaining_nodes.remove(nb1)
            remaining_nodes.remove(nb2)
            for node in remaining_nodes:
                dist = 0.5 * (d(nb1, node) + d(nb2, node) - d(nb1, nb2))
                distmat.setdist(newnode, node, dist)
            remaining_nodes.append(newnode)

            # (3) Update summed dists
            # For each node compute change from previous value and alter accordingly
            # Note: only nodes in "remaining_nodes" list are considered, even though "distmat" may contain more
            udist[newnode] = 0.0
            for node in remaining_nodes:

                # Add current value to new entry for "new"
                newdist = d(node, newnode)
                udist[newnode] += newdist

                # Update old entry for "node"
                diff = newdist - d(node, nb1) - d(node, nb2)
                udist[node] += diff

        # After loop. Set length of branch conecting final two nodes
        n1, n2 = remaining_nodes[0], remaining_nodes[1]
        dist = d(n1, n2)
        njtree.deroot()
        njtree.setlength(n1, n2, dist)

        return njtree

    #######################################################################################

    def upgma(self):
        """Computes UPGMA tree, returns Tree object"""

        # Depth of node = distance from leaf-level to node (i.e., depth of leaves = 0.0)
        depth = {leaf:0.0 for leaf in self.leaves}

        # Construct star-tree. This will be resolved node by node during algorithm
        upgmatree = Tree.from_leaves(self.leaves)
        rootnode = upgmatree.root

        # Construct copy of distmatrix: this will be altered during run (original will remain unchanged so more trees can be made)
        distmatcopy = copy.deepcopy(self.distmat)

        # Main loop: continue merging nearest nodes, on tree and in distmat, until only two nodes left
        while len(distmatcopy) > 1:
            (dist, (node1, node2)) = distmatcopy.nearest()

            # insert new node below two nodes to be merged, unless only two nodes remain (in which case: use root node)
            if len(distmatcopy) == 2:
                mergenode = rootnode
            else:
                mergenode = upgmatree.insert_node(rootnode, [node1, node2])

            depth[mergenode] = 0.5 * dist
            dist1 = depth[mergenode] - depth[node1]
            dist2 = depth[mergenode] - depth[node2]
            upgmatree.setlength(mergenode, node1, dist1)
            upgmatree.setlength(mergenode, node2, dist2)
            distmatcopy.merge_nodes(node1, node2, mergenode)

        return upgmatree

#############################################################################################
#############################################################################################


# # Placeholder: Insert test code here and run module in standalone mode
def main():
    pass

#######################################################################################


if __name__ == "__main__":
    main()

########################
#
##  Code for function profiling:
    # import cProfile
    # cProfile.run('main()', 'mainprofile')
    #
##    In separate window do:
##    >>> import pstats
##    >>> p = pstats.Stats('mainprofile')
##    >>> p.sort_stats('time').print_stats(10)
########################
#
# line_profiler:
# 1) Add @profile decorator just above def statement for function(s?) you would like to be profiled
# 2) Run these commands in shell:
#     kernprof -l treelib_devel.py
#     /usr/local/bin/python3 -m line_profiler treelib_devel.py.lprof
########################
#
# Memory profiler:
# 1) Add @profile decorator just above def statement for function(s?) you would like to be profiled
# 2) /usr/local/bin/python3 -m memory_profiler treelib_devel.py
#
# plot over time:
# mprof run treelib_devel.py
# mprof plot
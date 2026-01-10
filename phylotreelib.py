"""Classes and methods for analyzing, manipulating, and building phylogenetic trees"""
# Anders Gorm Pedersen
# Section for Bioinformatics, DTU Health Technology, Technical University of Denmark
# agpe@dtu.dk
import copy
import functools
import itertools
from itertools import (takewhile,repeat)
import math
import random
import re
import statistics
import sys
from io import StringIO
from operator import itemgetter
from collections import Counter
from collections import defaultdict
import numpy as np
#from line_profiler import profile

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

####################################################################################

def fast_treecount(filename, fileformat="nexus"):
    """Heuristic: count patterns ([;=\n] etc) to infer number of trees"""

    # Empirically: if ); is in file, then this == number of trees
    n_terminators = count_bytestring(filename, b");")
    if n_terminators > 0:
        return n_terminators

    # Count semicolon: if == 1 (possibly after nexus correction): 1 tree
    n_semicolons = count_bytestring(filename, b";")
    if n_semicolons == 1:
        return 1
    if fileformat == "nexus":
        n_other_semicolon_patterns  = count_bytestring(filename, b"egin taxa;")
        n_other_semicolon_patterns += count_bytestring(filename, b"egin trees;")
        n_other_semicolon_patterns += count_bytestring(filename, b"nd;")
        n_other_semicolon_patterns += count_bytestring(filename, b"imensions ntax")
        n_other_semicolon_patterns += count_bytestring(filename, b"axlabels")
        n_other_semicolon_patterns += count_bytestring(filename, b"ranslate")
        n_semicolons -= n_other_semicolon_patterns
    if n_semicolons == 1:
        return 1

    # If we got this far, and filetype is newick then bail out and use precise counting
    if fileformat == "newick":
        return count_trees_by_parsing(filename, fileformat)

    # Final attempt to infer ntrees for nexus files:
    # count "= (", and "tree "
    # Add the values that are not 0 to list, choose minimum as count
    # Should be robust to most variations, but should check at end of sumt...
    n_eqparen = count_bytestring(filename, b"= (")
    n_treestr = count_bytestring(filename, b"tree ")
    countlist = [n_semicolons, n_eqparen, n_treestr]
    notzero = [val for val in countlist if val>0]
    return min(notzero)

####################################################################################

def count_bytestring(filename, bytestring):
    """Fast counting of specific pattern. Bytestring argument must be given
    with b modifier (e.g., b');')"""

    # Modified from: https://stackoverflow.com/a/27517681/7836730
    with open(filename, 'rb') as f:
        bufsize = 1024*1024
        bufgen = takewhile(lambda x: x, (f.raw.read(bufsize) for _ in repeat(None)))

        prev_buf = b""
        count = 0

        for buf in bufgen:
            count += buf.count(bytestring)

            # For multi-byte patterns, consider overlaps between buffers
            if len(bytestring) > 1 and len(prev_buf) > 0:
                merged = prev_buf[-len(bytestring)+1:] + buf[:len(bytestring)-1]
                count += merged.count(bytestring)

            prev_buf = buf

    return count

####################################################################################

def count_trees_by_parsing(filename, args):
    # Open treefile. Discard (i.e., silently pass by) the requested number of trees
    if args.informat == "nexus":
        treefile = Nexustreefile(filename)
    else:
        treefile = Newicktreefile(filename)
    treecount = 0
    for tree in treefile:
        treecount += 1
    return treecount

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
        self.leaftup_interndict = {}
        self.unhashable_dict = {}

    def intern_leafset(self, leafset):
        return self.leafset_interndict.setdefault(leafset, leafset)

    def intern_clade(self, clade):
        return self.clade_interndict.setdefault(clade, clade)

    def intern_bipart(self, bipart):
        return self.bip_interndict.setdefault(bipart, bipart)

    def intern_leaftup(self, leaftup):
        return self.leaftup_interndict.setdefault(leaftup, leaftup)

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

    __slots__ = ("length", "label",
                 "SUMW", "mean", "M2", "n",
                 "length_var", "length_sd",
                 "posterior", "freq", "bipartition_cred", "rootcred",
                 "parent_height", "kid_height")

    # Used by __str__ when deciding which attributes to print (those with non-defaul values)
    _DEFAULTS = {
        "length": 0.0,
        "label": "",
        "SUMW": 0.0,
        "mean": 0.0,
        "M2": 0.0,
        "n": 0,
        "length_var": 0.0,
        "length_sd": 0.0,
        "posterior": None,
        "freq": None,
        "bipartition_cred": None,
        "rootcred": None,
        "parent_height": None,
        "kid_height": None,
    }

    def __init__(self, length=0.0, label=""):
        self.length = length
        self.label = label
        self.SUMW = 0.0
        self.mean = 0.0
        self.M2 = 0.0
        self.n = 0
        self.length_var = 0.0
        self.length_sd = 0.0
        self.posterior = None
        self.freq = None
        self.bipartition_cred = None
        self.rootcred = None
        self.parent_height = None
        self.kid_height = None

    ###################################################################################################
    # The code below is used by Tree.__str__ when choosing which attr to print)

    @classmethod
    def printable_attrs(cls):
        """Attrs eligible to appear as extra columns (excludes 'length' because it's printed as branchlen)."""
        return tuple(a for a in cls.__slots__ if a != "length")

    @classmethod
    def default_for(cls, attr: str):
        return cls._DEFAULTS.get(attr, None)

    @classmethod
    def is_interesting(cls, attr: str, value) -> bool:
        """Decide whether a value should cause the column to be included / printed."""
        default = cls.default_for(attr)
        return value != default

    def interesting_attrs(self):
        """Yield attrs that are interesting on *this* branch (non-default)."""
        for attr in self.printable_attrs():
            v = getattr(self, attr)
            if self.is_interesting(attr, v):
                yield attr

    def format_attr(self, attr: str) -> str:
        """String to show in table cell (boring values become '')."""
        v = getattr(self, attr)
        if not self.is_interesting(attr, v):
            return ""
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)

    ###################################################################################################


    def __str__(self):
        return f"length: {self.length}\n"

    def __repr__(self):
        return self.__str__()

    ###############################################################################################

    def copy(self, copylengths=True, copyattr=True):
        """Returns copy of Branchstruct object, with selected attributes included
        copylength: copy .length attribute
        copyattr: copy remaining attributes"""

        # Python note: only works as deepcopy if all attributes are immutable. Rewrite!
        obj = Branchstruct()
        if copylengths:
            obj.length = self.length
        if copyattr:
            for name in self.__slots__[1:]:
                setattr(obj, name, getattr(self, name))
        return obj

    ###############################################################################################

    def merge(self, other, check_compat=False):
        """Merges two Branchstructs and returns new Branchstruct, e.g. for collapsing root branches.
        lengths are summed, all other attributes are ignored (no easy way to combine)"""

        # Python note: non-length attributes difficult to manage, since they may differ on two
        # branches, and there is no way of combining them (e.g. clade cred-derived support values)
        # Consider whether this will ever impact callers

        obj = Branchstruct(length = self.length + other.length)
        return obj

###################################################################################################
###################################################################################################

class Nodestruct:
    """Class that emulates a struct. Keeps node-related info"""

    __slots__ = ("depth", "nleaves", "subcladepairs", "best_pair",
                 "SUMW", "mean", "M2", "n",
                 "clade_cred", "posterior", "freq",
                 "depth_var", "depth_sd",
                 "clade_score")

    def __init__(self, depth=0.0, nleaves=None):
        self.depth = depth
        self.nleaves = nleaves
        self.subcladepairs = set()
        self.best_pair = None
        self.SUMW = 0.0
        self.mean = 0.0
        self.n = 0
        self.clade_cred = None
        self.posterior = None
        self.freq = None
        self.depth_var = None
        self.depth_sd = None
        self.clade_score = None

    def __str__(self):
        return f"depth: {self.depth}\n"

    def __repr__(self):
        return self.__str__()

    ###############################################################################################

    def add_subcladepair(self, subclade1, subclade2):
        self.subcladepairs.add(frozenset([subclade1, subclade2]))

    ###############################################################################################

    def copy(self):
        """Returns copy of Nodestruct object, with all attributes included"""

        # Python note: will cause problems if any values are not immutable
        obj = Nodestruct()
        for name in self.__slots__:
            setattr(obj, name, getattr(self, name))
        return obj

###################################################################################################
###################################################################################################

class Topostruct:
    """Class that emulates a struct. Keeps topology-related info"""

    __slots__ = ["weight", "tree", "posterior"]

    # Python note: perhaps replace with dataclass, available since python 3.7
    pass

###################################################################################################
###################################################################################################

class Bipartition:
    """Class that represents a bipartition: the split of leaves in two sets that corresponds
    to a branch in a tree. Class has a fast method for hashing and therefore useful when
    comparing Bipartitions"""

    __slots__ = ("_mask", "_sorted_leaf_tup")

    # Class-level cache
    # dict of {_tipset_id: object_cache}
    # object_cache: {canonical_mask: Bipartition}
    _class_cache = {}

    # Called only from other constructors, not directly
    # mask: int interpreted as bitset of leaf indices
    def __init__(self, mask, sorted_leaf_tup):
        self._mask = mask
        self._sorted_leaf_tup = sorted_leaf_tup

    @classmethod
    def from_halfmask(cls, halfmask, alltips_mask, object_cache, sorted_leaf_tup):
        """Bipartition constructor that relies on caller knowing leaf universe and providing
        required information:

        halfmask: int representing bitset of indices for half of bipartition
        object_cache: class-level cache for bipartition objects belonging to this leaf-universe
        alltips_mask, sorted_leaf_tup: for this leaf universe; cached and interned on Tree class

        returns Bipartition object
        Achieves interning + caching to avoid repeated object construction"""

        # Represent bipartition using just one half of leaves
        # Choose bitset corresponding to smaller int (canonical_mask)
        canonical_mask = min(halfmask, alltips_mask ^ halfmask)
        obj = object_cache.get(canonical_mask)
        if obj is None:
            obj = cls(canonical_mask, sorted_leaf_tup)
            object_cache[canonical_mask] = obj
        return obj

    @classmethod
    def from_halfmask_unknown_leafuniverse(cls, halfmask, tree):
        """halfmask: int representing bitset of indices for half of bipartition
        tree: Tree object that this Bipartition belongs to (for access to cached attributes)
        returns Bipartition object
        Achieves interning + caching to avoid repeated object construction"""

        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = tree.cached_attributes

        object_cache = cls._class_cache.get(frozenset_leaves)
        if object_cache is None:
            object_cache = {}
            cls._class_cache[frozenset_leaves] = object_cache

        canonical_mask = min(halfmask, halfmask ^ alltips_mask)
        obj = object_cache.get(canonical_mask)
        if obj is None:
            obj = cls(canonical_mask, sorted_leaf_tup)
            object_cache[canonical_mask] = obj
        return obj

    @classmethod
    def from_leafset(cls, leafset, tree):
        # optional compatibility constructor if some callers still pass leaf names
        halfmask = 0
        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = tree.cached_attributes
        for leaf in leafset:
            #halfmask |= (1 << leaf2idx[leaf])
            halfmask += leaf2mask[leaf]
        return cls.from_halfmask_unknown_leafuniverse(halfmask, tree)

    @staticmethod
    def _canonical_half(halfmask, alltips_mask, ntips):
        """Return canonical representation of a split (smaller side; deterministic tie-break).
        masks are ints representing bitsets for leaf indices"""
        k = halfmask.bit_count()
        if k > ntips - k:
            halfmask = alltips_mask ^ halfmask
            k = ntips - k
        elif k == ntips - k:
            other = alltips_mask ^ halfmask
            if other < halfmask:
                halfmask = other
        return halfmask

    def __hash__(self):
        # Python note: this strategy (using only the halfmask as hash) depends on never comparing
        # Bipartitions across leaf-universes. In class-level cache this is ensured by having
        # object_cache be indexed by tipset_id
        return self._mask

    def __eq__(self, other):
        # Python note: see note under __hash__
        if self is other:
            return True
        return self._mask == other._mask

    def get_bipartitions(self):
        """Return (side1, side2) as frozensets of leaf names."""
        sorted_leaf_tup = self._sorted_leaf_tup
        ntips = len(sorted_leaf_tup)
        side1 = frozenset(sorted_leaf_tup[i] for i in range(ntips) if (self._mask >> i) & 1)
        side2 = frozenset(sorted_leaf_tup) - side1
        return side1, side2
        
    def get_canonical_clade(self):
        """Returns a frozenset of leaves included in canonical mask (= side1 from get_bipartitions)"""
        side1,side2 = self.get_bipartitions()
        return side1
        
    def __iter__(self):
        return iter(self.get_bipartitions())

    def __str__(self):
        bip1,bip2 = self.get_bipartitions()
        if len(bip1) < len(bip2):
            return f"leaf set 1:\n{str(bip1)}\nleaf set 2:\n{str(bip2)}"
        else:
            return f"leaf set 1:\n{str(bip2)}\nleaf set 2:\n{str(bip1)}"

    def __repr__(self):
        return self.__str__()

###################################################################################################
###################################################################################################

class Clade:
    """Class that represents a clade: the set of leaves descended from an internal node"""

    __slots__ = ("_mask", "_sorted_leaf_tup")

    # Class-level cache
    # dict of {_tipset_id: object_cache}
    # object_cache: {mask: Clade}
    _class_cache = {}

    # Called only from other constructors, not directly
    # mask: int interpreted as bitset of leaf indices in clade
    def __init__(self, mask, sorted_leaf_tup):
        self._mask = mask
        self._sorted_leaf_tup = sorted_leaf_tup

    @classmethod
    def from_mask(cls, mask, object_cache, sorted_leaf_tup):
        """Clade constructor that relies on caller knowing leaf universe and providing
        required information:

        mask: int representing bitset of indices for leaves in Clade
        object_cache: class-level cache for clade objects belonging to this leaf-universe
        sorted_leaf_tup: for this leaf universe; cached and interned on Tree class

        returns Clade object
        Achieves interning + caching to avoid repeated object construction"""

        cached_clade = object_cache.get(mask)
        if cached_clade is None:
            cached_clade = cls(mask, sorted_leaf_tup)
            object_cache[mask] = cached_clade
        return cached_clade

    @classmethod
    def from_leafset(cls, leafset, tree):
        """Convenience constructor for tests and interactive use where clade is given as
        leafset, and this method has to handle object_cache retrieval and construction."""

        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = tree.cached_attributes

        object_cache = cls._class_cache.get(frozenset_leaves)
        if object_cache is None:
            object_cache = {}
            cls._class_cache[frozenset_leaves] = object_cache

        mask = 0
        for leaf in leafset:
            mask += leaf2mask[leaf]

        cached_clade = cls.from_mask(mask, object_cache, sorted_leaf_tup)

        return cached_clade

    def __hash__(self):
        # Python note: this strategy (using only the mask as hash) depends on never comparing
        # Clades across leaf-universes. In class-level cache this is ensured by having
        # object_cache be indexed by tipset_id
        return self._mask

    def __eq__(self, other):
        if self is other:
            return True
        return self._mask == other._mask

    def __len__(self):
        return len(self._sorted_leaf_tup)

    def get_clade(self):
        """Return leaf names in this clade."""
        sorted_leaf_tup = self._sorted_leaf_tup
        ntips = len(sorted_leaf_tup)
        leafnames = frozenset(sorted_leaf_tup[i] for i in range(ntips) if (self._mask >> i) & 1)
        return leafnames

    # Python note: this allows unpacking as if the class was a tuple: c1 = myclade
    # 1-tuple; comma is important
    def __iter__(self):
        return iter((self.get_clade(),))

    def __str__(self):
        return f"\n{self.get_clade()}\n"

    def __repr__(self):
        return self.__str__()

###################################################################################################
###################################################################################################

class TreeError(Exception):
    pass

###################################################################################################
###################################################################################################

class NewickStringParser:
    """Class creating parser for specific Newick tree string. Used in Tree.from_string"""

    def __init__(self, transdict=None, label_attr_name="label", label_type=str):
        """Initialise parser object.
        label_attr_name: name of label attribute on Branchstructs
        label_type: type conversion to perform on input label string (e.g. float)

        To create derived class:
          (1) Call super __init__ in __init__ of derived class
          (2) Override set_token_delimiters() (return string of single-char delimiters)
          (3) Override update_dispatch() to modify self.dispatch dict"""

        self.delimiters = self.set_token_delimiters()
        self.delimset = set(self.delimiters)
        self.regex = self.create_compiled_regex(self.delimiters)

        # If transdict was supplied: use corresponding handler when adding leaves
        # Also store transdict in parser (at class level! this is so staticmethods can access)
        if transdict:
            self.transdict = transdict
            self._handle_add_leaf = self._handle_add_leaf_with_transdict

        # Control how labels are interpreted and what the Branchstruct attribute will be named
        self.label_attr_name = label_attr_name
        self.label_type = label_type

        # Dispatch dictionary:
        # Key = current state and token-type
        # Value = tuple of
        # (1) Function that should be run to further build tree data structure (argument=token-value)
        # (2) New state (the one we move to after this step)
        self.dispatch = {
            "TREE_START":       {   "(":        (self._handle_add_root_intnode,     "INTNODE_START")    },
            "INTNODE_START":    {   "(":        (self._handle_add_intnode,          "INTNODE_START"),
                                    "NUM_NAME": (self._handle_add_leaf,             "LEAF")             },
            "LEAF":             {   ":":        (self._handle_transition_brlen,     "EXPECTING_BRLEN"),
                                    ",":        (self._handle_transition_child,     "EXPECTING_CHILD"),
                                    ")":        (self._handle_intnode_end,          "INTNODE_END")      },
            "EXPECTING_BRLEN":  {   "NUM_NAME": (self._handle_add_brlen,            "BRLEN")            },
            "BRLEN":            {   ",":        (self._handle_transition_child,     "EXPECTING_CHILD"),
                                    ")":        (self._handle_intnode_end,          "INTNODE_END")      },
            "EXPECTING_CHILD":  {   "(":        (self._handle_add_intnode,          "INTNODE_START"),
                                    "NUM_NAME": (self._handle_add_leaf,             "LEAF")             },
            "INTNODE_END":      {   ")":        (self._handle_intnode_end,          "INTNODE_END"),
                                    ",":        (self._handle_transition_child,     "EXPECTING_CHILD"),
                                    ":":        (self._handle_transition_brlen,     "EXPECTING_BRLEN"),
                                    "NUM_NAME": (self._handle_label,                "LABEL"),
                                    ";":        (self._handle_transition_tree_end,  "TREE_START")       },
            "LABEL":            {   ")":        (self._handle_intnode_end,          "INTNODE_END"),
                                    ",":        (self._handle_transition_child,     "EXPECTING_CHILD"),
                                    ":":        (self._handle_transition_brlen,     "EXPECTING_BRLEN")  }
        }

        self.update_dispatch()

    ###############################################################################################

    def set_token_delimiters(self):
        """Override to modify the set of delimiters used when parsing.
        All delimiters must be single-character (handle multi-char delims in parsing code)
        Specify single string of delimiters that is returned"""

        return ",:;()"          # Override this line in derived classes

    ###############################################################################################

    def update_dispatch(self):
        """Override to modify dispatch dictionary (extra edges and nodes in state diagram)"""
        pass        # Replace this line with code modifying dispath dict in derived classes
                    # For instance: self.dispatch["INTNODE_END"]["["] = self._handle_enter_comment

    ###############################################################################################

    def create_compiled_regex(self, delimiters):
        """Creates compiled regex for faster parsing"""

        # Final pattern will be of the form "([abcde])": a capturing parenthesis surrounding
        # a bracketed set of single-char delimiters (here a, b, c, d, e)
        # Metachars are escaped using re.escape
        pattern_list = ["(["]
        for delim in delimiters:
            pattern_list.append(f"{re.escape(delim)}")
        pattern_list.append("])")
        pattern = "".join(pattern_list)
        return re.compile(pattern)

    ###############################################################################################

    def parse(self, treeobj, treestring):
        # Construct Tree object that is filled out while parsing
        # Tree is represented as a dictionary of dictionaries.
        # The keys in the top dictionary are the internal nodes which are numbered consecutively.
        # Each intnode key has an associated value that is itself a dictionary listing the children:
        # Keys are child nodes, values are Branchstructs which contain two attributes:
        #   (1) "length" (float, parsed from numbers after colons in input newick string)
        #   (2) label_attr_name (by default "label", but is specified at runtime) from before :
        #       input string for this attribute is converted using label_type() (e.g. float or str)
        # Leafs are identified by a string instead of a number

        # NOTE: interprets non-leaf labels as belonging to an internal branch (not to
        # an internal node). The label is attached to the same branch as the branch length
        self.treeobj = treeobj
        self.treeobj.root = 0

        # These variables are only used while parsing, and should be in parserobj (not treeobj)
        # Values are reset each time new treestring is parsed
        self.node_stack = []
        self.ntips = 0
        self.nodeno = None

        treestring = "".join(treestring.split()) # Hack to remove whitespace from treestring
        dispatch = self.dispatch
        delimset = self.delimset
        state = "TREE_START"
        tree_parts_list = self.regex.split(treestring)
        for token_value in tree_parts_list:
            if token_value:
                if token_value in delimset:
                    token_type = token_value
                else:
                    token_type = "NUM_NAME"
                try:
                    handler, state = dispatch[state][token_type]
                except KeyError:
                    self._handle_parse_error(state, token_value, token_type, treestring, tree_parts_list)
                handler(token_value)

        self.sanitychecks(self.treeobj, treestring)

        return self.treeobj

    ###############################################################################################

    def sanitychecks(self, treeobj, treestring):
        if self.node_stack:
            msg = f"Nodes remaining on stack after parsing treestring: {self.node_stack}"
            raise TreeError(msg)
        if self.ntips > len(treeobj.leaves):
            raise TreeError(f"Duplicated leafnames in treestring:\n{treestring}")

    ###############################################################################################

    def _handle_parse_error(self, state, token_value, token_type, treestring, tree_parts_list):
        # If unexpected token-type was encountered: first check if parentheses are balanced:
        if treestring.count("(") != treestring.count(")"):
            msg = "Imbalance in tree-string: different number of left- and right-parentheses\n"
            msg += f"Left: ({treestring.count('(')}  Right: {treestring.count(')')})"
            raise TreeError(msg)
        else:
            # If that was not the problem: report current parser-state, token-type, and token-value
            msg = ("Parsing error: unexpected token-type for this state:\n"
                   f"Parser state:     {state}\n"
                   f"Token_type:       {token_type}\n"
                   f"Token-value:      {token_value}\n"
                   f"Tree-string:      {treestring}\n"
                   f"Tree-parts list:  {tree_parts_list}\n")
                   
            raise TreeError(msg)

    ###############################################################################################

    def _handle_add_root_intnode(self, token_value):
        self.nodeno = 0
        self.treeobj.child_dict[self.nodeno] = {}
        self.node_stack.append(self.nodeno)

    ###############################################################################################

    def _handle_add_intnode(self, token_value):
        self.nodeno += 1
        self.treeobj.child_dict[self.nodeno] = {}
        parent = self.node_stack[-1]
        self.treeobj.child_dict[parent][self.nodeno] = Branchstruct()
        self.node_stack.append(self.nodeno)

    ###############################################################################################

    def _handle_add_leaf(self, name):
        child = sys.intern(name)
        parent = self.node_stack[-1]
        self.treeobj.child_dict[parent][child] = Branchstruct()
        self.node_stack.append(child)
        self.treeobj.leaves.add(child)
        self.ntips += 1

    ###############################################################################################

    def _handle_add_leaf_with_transdict(self, name):
        child = sys.intern(self.transdict[name])
        parent = self.node_stack[-1]
        treeobj = self.treeobj
        treeobj.child_dict[parent][child] = Branchstruct()
        self.node_stack.append(child)
        treeobj.leaves.add(child)

    ###############################################################################################

    def _handle_transition_child(self, token_value):
        self.node_stack.pop()

    ###############################################################################################

    def _handle_transition_brlen(self, token_value):
        pass

    ###############################################################################################

    def _handle_add_brlen(self, brlen_string):
        try:
            brlen = float(brlen_string)
            child = self.node_stack[-1]
            parent = self.node_stack[-2]
            self.treeobj.child_dict[parent][child].length = brlen
        except ValueError:
            raise TreeError(f"Expected branch length: {brlen_string}")

    ###############################################################################################

    def _handle_intnode_end(self, token_value):
        self.node_stack.pop()

    ###############################################################################################

    def _handle_label(self, label):
        child = self.node_stack[-1]
        parent = self.node_stack[-2]
        branchstruct = self.treeobj.child_dict[parent][child]
        setattr(branchstruct, self.label_attr_name, self.label_type(label))

    ###############################################################################################

    def _handle_transition_tree_end(self, token_value):
        self.node_stack.pop()

###################################################################################################
###################################################################################################

class Tree:
    """Class representing basic phylogenetic tree object."""

    # Tree is represented as a dictionary of dictionaries. The keys in the top dictionary
    # are the internal nodes, which are numbered consecutively. Each key has an associated
    # value that is itself a dictionary listing the children: keys are child nodes, values are
    # Branchstructs containing "length" (float) and "label" (str) fields.
    # Leafs are identified by a string instead of a number

    # Implementation note: Tree objects can be constructed from several different kinds of things:
    # including Newick tree strings, Bipartition lists, and a list of leaves.
    # The Tree class therefore has several alternate constructors implemented as classmethods
    # The main constructor "__init__" is therefore mostly empty

    # Class-level cache;
    # dict of: {frozenset_leaves: (frozenset_leaves, sorted_leaf_tup, alltips_mask, ntips) }
    # Key is tipset (frozenset_leaves), which functions as a universe-id (in case there are
    # different types of trees in one Python session)
    # Cache achieves interning (one copy per treetype of the variables in the tuple)
    # and also avoiding repeated construction of the variables
    _class_cache = {}

    def __init__(self):
        self.child_dict = {}
        self._nodedict = None
        self.leaves = set()
        self._intnodes = None
        self._nodes = None
        self._parent_dict = None
        self._remotechildren_dict = None
        self._remotechildren_indices_dict = None
        self._remotechildren_mask_dict = None
        self._frozenset_leaves = None
        self._leaf2index = None
        self.dist_dict = None
        self._pathdist_dict = None
        self._pathdist_dict_unroot = None
        self._pathdist_as_ndarray = None
        self._pathdist_as_ndarray_unroot = None
        self.path_dict = None
        self.interner = None
        self._sorted_intnodes_deep = None
        self._sorted_intnodes_shallow = None
        self._rootdist = None
        self._nodedepthdict = None
        self._topology_bipart = None
        self._topology_clade = None
        self._cached_attributes = None

    ###############################################################################################

    def clear_caches(self, preserve=[]):
        """Clears set of computed attributes. Use when tree structure is changed, and these
        are no longer reliable. Does not clear attributes listed in preserve"""

        # structure/topology derived
        for attr in {
            "_parent_dict", "_remotechildren_dict", "_remotechildren_indices_dict",
            "_remotechildren_mask_dict",
            "_cached_attributes", "_leaf2index",
            "path_dict", "_intnodes", "_nodes",
            "_sorted_intnodes_deep", "_sorted_intnodes_shallow",
            "_topology_bipart", "_topology_clade",
        } - set(preserve):
            setattr(self, attr, None)

        # and also clear length-derived
        self.clear_length_caches()

    ###############################################################################################

    def clear_length_caches(self):
        """Clear only length-related caches - leave topology related caches alone"""

        for attr in {
            "dist_dict", "_rootdist", "_nodedepthdict", "_pathdist_dict",
            "_pathdist_dict_unroot", "_pathdist_as_ndarray", "_pathdist_as_ndarray_unroot",
        }:
            setattr(self, attr, None)
        self.nodedist.cache_clear()

    ###############################################################################################

    @classmethod
    def from_string(cls, orig_treestring, transdict=None, format="newick",
                    label_attr_name="label", label_type=str):
        """Constructor: Tree object from tree-string in specified format"""

        if format.lower() == "newick":
            parser_obj = NewickStringParser(transdict, label_attr_name, label_type)
        treeobj = cls._from_string_private(parser_obj, orig_treestring)
        del parser_obj
        return treeobj

    ###############################################################################################

    @classmethod
    def _from_string_private(cls, parser_obj, orig_treestring):
        """Constructor: Tree object from tree-string.
        This version is meant to be used by e.g. treefile objects that will use same parser
        for a lot of tree-strings, and can then instantiate parser once and reuse with this method.
        Instantiated parser object has to have same format as treestring (for instance
        Newick or Nexus with metacomments)"""

        # All action is in parser_obj
        treeobj = cls()
        treeobj = parser_obj.parse(treeobj, orig_treestring)
        return treeobj

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
    def from_biplist(cls, biplist):
        """Constructor: Tree object from bipartition list"""

        # Input is a bipartitionlist (actually a dictionary of bipartition:Branchstruct pairs):
        # Names of leaves on one side of a branch are represented as an immutable set of leaves
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree is represented as a dictionary of bipartition:Branchstruct pairs

        # Extract set of leaves
        part1, part2 = next(iter(biplist))  # First key=set of two leaf name sets
        leaves = part1 | part2          # Concatenate them to get all leafnames

        # Start by building star-tree, this is resolved branch-by-branch later on
        obj = Tree.from_leaves(leaves)

        # Iterate over all bipartitions, for each: add extra branch and/or update Branchstruct
        for (bip1, bip2), branchstruct in biplist.items():

            # If bipartition represents external branch: update relevant Branchstruct
            if len(bip1) == 1 or len(bip2) == 1:
                if len(bip1) == 1:
                    (leaf, ) = bip1      # A one-member tuple used for value unpacking (pop)
                else:
                    (leaf, ) = bip2

                # Find child_edges dict containing leaf, and update branchstruct
                for child_edges in obj.child_dict.values():
                    if leaf in child_edges:
                        child_edges[leaf] = branchstruct
                        break

            # If bipartition represents internal branch:
            # insert new node and attach branchstruct there
            else:
                # Make bip1 the smaller part (for MRCA efficiency)
                if len(bip1) > len(bip2):
                    bip1,bip2 = bip2,bip1

                # Determine which group of leaves to move
                # Note: one part of bipartition will necessarily have root as its MRCA
                # since the members are present on both sides of root. It is the other part of
                # the bipartition (where all members are on same side of root) that should be moved
                # For a star-tree the resolution will be random (both have root as their MRCA)
                mrca1 = obj.find_mrca(bip1)
                if mrca1 != obj.root:
                    insertpoint = mrca1
                    active_bip = bip1
                else:
                    mrca2 = obj.find_mrca(bip2)
                    insertpoint = mrca2
                    active_bip = bip2

                # Determine which of insertpoint's children to move, namely all children
                # that are either in the active bipartition OR children whose descendants are
                movelist = []
                for child in obj.children(insertpoint):
                    if child in active_bip:
                        movelist.append(child)
                    elif obj.remotechildren_dict[child] <= active_bip:
                        movelist.append(child)

                # Add new internal node and move relevant children to this.
                # Attach Branchstruct to newly created branch
                # Note: this method calls clear_caches
                obj.insert_node(insertpoint, movelist, branchstruct)

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

        # Extract set of leaves
        part1, part2 = next(iter(topology)) # First item=set of two leaf name sets
        leaves = part1 | part2          # Concatenate them to get all leafnames

        # Start by building star-tree, this is resolved branch-by-branch later on
        obj = Tree.from_leaves(leaves)

        # Iterate over all bipartitions, for each: add extra branch and/or update Branchstruct
        for bip1, bip2 in topology:

            # If bipartition represents internal branch: add branch to tree
            if len(bip1) > 1 and len(bip2) > 1:

                # Make the smaller side bip1 for MRCA efficiency
                if len(bip1) > len(bip2):
                    bip1,bip2 = bip2,bip1

                # Determine which group of leaves to move
                # One part of bipartition will necessarily have root as its MRCA
                # since the members are present on both sides of root. It is the other part of
                # the bipartition (where all members are on same side of root) that should be moved
                # For a star-tree the resolution will be random (both have root as their MRCA)
                mrca1 = obj.find_mrca(bip1)
                if mrca1 != obj.root:
                    insertpoint = mrca1
                    active_bip = bip1
                else:
                    mrca2 = obj.find_mrca(bip2)
                    insertpoint = mrca2
                    active_bip = bip2

                # Determine which of insertpoint's children to move, namely all children
                # that are either in the active bipartition OR children whose descendants are
                movelist = []
                for child in obj.children(insertpoint):
                    if child in active_bip:
                        movelist.append(child)
                    elif obj.remotechildren_dict[child] <= active_bip:
                        movelist.append(child)

                # Add new internal node and move relevant children to this.
                # Attach empty Branchstruct to newly created branch (no blen info in topology)
                obj.insert_node(insertpoint, movelist, Branchstruct())

        return obj

    ###############################################################################################

    @classmethod
    def from_cladedict(cls, cladedict):

        clade = next(iter(cladedict))
        obj = cls.from_leaves(clade._sorted_leaf_tup)
        all_leaves = obj.frozenset_leaves

        for clade, nodestruct in cladedict.items():
            clade_leaves = clade.get_clade()
            if len(clade_leaves) == 1:
                leaf = next(iter(clade_leaves))
                obj.nodedict[leaf] = nodestruct
            elif clade_leaves == all_leaves:
                obj.nodedict[obj.root] = nodestruct
            else:
                mrca = obj.find_mrca(clade_leaves)
                movelist = []
                for child in obj.children(mrca):
                    if (child in clade_leaves) or (obj.remotechildren_dict[child] <= clade_leaves):
                        movelist.append(child)
                newnode = obj.insert_node(mrca, movelist, Branchstruct())
                obj.nodedict[newnode] = nodestruct

                # Reset obj caches which are now obsolete
                obj.clear_caches()

        return obj

    ###############################################################################################

    @classmethod
    def from_branchinfo(cls, parentlist, childlist, lenlist=None, **attrlists):
        """Constructor: Tree object from information about all branches in tree

        Information about one branch is conceptually given as:
            parentnodeID, childnodeID, [length], [other named attributes...]

        The function takes as input 2 or more separate lists containing:
            IDs of parents (internal nodes, integer or string)
            IDs of children (internal or leaf nodes, integer or string)
            Length of branches (optional)
            Any number of additional attributes given as <keyword=list>

        The lists are assumed to have same length and be in same order (so index n in
        each list corresponds to same branch).

        Note: most IDs appear multiple times in lists
        Note 2: can be used as workaround so user can specify IDs for internal nodes"""

        nbranches = len(parentlist)
        if len(childlist) != nbranches:
            raise TreeError(f"List 'childlist' does not have same length as parentlist: {len(childlist)} != {len(parentlist)}")
        if lenlist is None:
            lenlist = [0.0]*nbranches
        elif len(lenlist) != nbranches:
            raise TreeError(f"List 'childlist' does not have same length as parentlist: {len(lenlist)} != {len(parentlist)}")
        for attrname, attrlist in attrlists.items():
            if len(attrlist) != nbranches:
                raise TreeError(f"List '{attrname}' does not have same length as parentlist: {len(attrlist)} != {len(parentlist)}")
        obj = cls()                    # Ensures class will be correct also for subclasses of Tree

        for i in range(nbranches):
            parent = parentlist[i]     # Perhaps check types are OK?
            child = childlist[i]
            blen = lenlist[i]
            if attrlists:
                attrdict = {name: value[i] for name,value in attrlists.items()}
            else:
                attrdict = dict()
            branch = Branchstruct(blen, **attrdict)
            if parent in obj.child_dict:
                obj.child_dict[parent][child] = branch
            else:
                obj.child_dict[parent] = { child:branch }

        # Leaves are the childnodes that are not in parentlist
        obj.leaves = set(childlist) - set(parentlist)

        # Root node is the parent node that is not also in childlist
        diffset = set(parentlist) - set(childlist)
        obj.root = diffset.pop()

        # Sanity check: are there any non-root, internal nodes that have no parent?
        # This would mean that sub-tree is not linked to rest of tree structure (missing branches)
        nonroot_intnodes = obj.intnodes - set([obj.root])
        orphans = nonroot_intnodes - set(childlist)
        if orphans:
            msg = (f"Missing branch-information: these (non-root) internal nodes have no parent: "
                   f"{list(orphans)}")
            raise TreeError(msg)

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
        """Prints table of parent-child relationships including branch lengths and other attributes."""

        # Parameters
        fixed_columns = ["parent", "child", "branchlen"]
        padding = 2  # Minimum padding spaces before and after each value
        max_col_width = 200
        max_total_width = 120  # Maximum width for the entire table

        # Collect data rows and attribute names
        data_rows = []
        attr_set = set()
        for parent in self.sorted_intnodes():
            kidset = self.children(parent)
            sorted_kids = sorted(kidset, key=lambda x: (isinstance(x, str), x))
            for kid in sorted_kids:
                parentstr = str(parent)
                kidstr = str(kid)
                branch = self.child_dict[parent][kid]
                dist = "{num:.6g}".format(num=branch.length)
                data_rows.append([parentstr, kidstr, dist, branch])

                # Slots-based: include only attrs that are non-default on at least one branch
                attr_set.update(branch.interesting_attrs())

        # Sort and initialize attribute names
        attr_list = sorted(attr_set)

        # Initialize headers and maxwidths for fixed columns
        headers = fixed_columns.copy()
        maxwidths = [max(len(col), 6) for col in fixed_columns]  # Python note: superseded by code below?

        # Compute max widths for fixed columns
        for i, col in enumerate(fixed_columns):
            maxw = len(col)
            for data_row in data_rows:
                value = data_row[i]
                maxw = min(max(len(value), maxw), max_col_width)
            maxwidths[i] = maxw

        # Compute max widths for attributes
        maxwidth_dict = {}
        for attr in attr_list:
            maxw = len(attr)
            for data_row in data_rows:
                branch = data_row[3]
                value_str = branch.format_attr(attr)
                maxw = min(max(len(value_str), maxw), max_col_width)
            maxwidth_dict[attr] = maxw

        # Include attributes until max_total_width is reached
        totwidth = sum(maxwidths) + len(headers) * (padding * 2 + 1) + 1
        included_attrs = []
        for attr in attr_list:
            potential_totwidth = totwidth + maxwidth_dict[attr] + (padding * 2 + 1)
            if potential_totwidth <= max_total_width:
                headers.append(attr)
                maxwidths.append(maxwidth_dict[attr])
                included_attrs.append(attr)
                totwidth = potential_totwidth
            else:
                break

        # Adjust totwidth for the table border
        totwidth = sum(maxwidths) + len(maxwidths) * (padding * 2 + 1) + 1

        # Build table rows
        table = [headers]
        for parentstr, kidstr, dist, branch in data_rows:
            row = [parentstr, kidstr, dist]
            for attr in included_attrs:
                value = getattr(branch, attr, "")
                value_str = branch.format_attr(attr)[:max_col_width]
                row.append(value_str)
            table.append(row)

        # Build the formatted table string
        border_line = "+" + "-" * (totwidth - 2) + "+\n"
        tabstring = border_line

        # Add header row
        for i, header in enumerate(headers):
            col_width = maxwidths[i] + padding * 2
            tabstring += "|" + (padding * " ") + header.center(col_width - 2*padding) + (padding * " ")
        tabstring += "|\n"
        tabstring += border_line

        # Add data rows
        for row in table[1:]:
            for i, value in enumerate(row):
                col_width = maxwidths[i] + padding * 2
                tabstring += "|" + (padding * " ") + str(value).rjust(col_width - 2 * padding) + (padding * " ")
            tabstring += "|\n"
        tabstring += border_line

        # Add list of leaves
        sorted_leaflist = sorted([str(leaf) for leaf in self.leaves])
        leaf_width = max(len(leaf) for leaf in sorted_leaflist)
        leaf_width = min(leaf_width, max_total_width)
        tabstring += f"\n{len(self.leaves)} Leaves:\n"
        tabstring += "-" * leaf_width + "\n"
        for leaf in sorted_leaflist:
            leaf_str = leaf[:leaf_width]
            tabstring += f"{leaf_str}\n"

        return tabstring

    ###############################################################################################

    def __eq__(self, other):
        """Implements unrooted equality testing for Tree objects"""

        return self.equals(other, rooted=False)

    ###############################################################################################

    def equals(self, other, rooted=False, blenprecision=0.005):
        """General method to test rooted or unrooted equality of Tree objects.
        Two trees are identical if they have the same leaves, the same topology
        and the same branchlengths (within requested relative precision).
        Branch labels are ignored.

        rooted=False: ignore rooting, compare .topology_bipart
        rooting=True: take rooting into account, compare .topology_clade
        blenprecision: floating point comparison of relative difference.
                       blens are identical if abs(len1 - len2) / len1) < blenprecision
                       Default precision 0.005 chosen based on empirical comparison
                       between own and PHYLIP trees (...)
        """

        if self.leaves != other.leaves:
            return False

        if rooted:
            topo_self = self.topology_clade
            topo_other = other.topology_clade
        else:
            topo_self = self.topology_bipart
            topo_other = other.topology_bipart
        if topo_self != topo_other:
            return False

        bipself = self.bipdict()
        bipother = other.bipdict()
        for bipart in bipself:
            len1 = bipself[bipart].length
            len2 = bipother[bipart].length
            if len1 != 0 and len2 != 0:
                if (abs(len1 - len2) / len1) > blenprecision:   # Floating point comparison of relative diff
                    return False
            if (len1 == 0 and len2 > 0) or (len1 > 0 and len2 == 0): # Will this ever happen?
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
    def cached_attributes(self):
        """Return tuple of cached, interned
                  (frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask)
        If not present on instance: Check class cache.
        If not present on class: build, cache on class and instance, and return"""

        if self._cached_attributes is None:
            frozenset_leaves = frozenset(self.leaves)
            class_cached_attr = type(self)._class_cache.get(frozenset_leaves)
            if class_cached_attr is None:
                sorted_leaf_tup = sorted(frozenset_leaves)
                leaf2index = {leaf:i for i,leaf in enumerate(sorted_leaf_tup)}
                leaf2mask = {leaf:(1<<i) for leaf,i in leaf2index.items()}
                ntips = len(sorted_leaf_tup)
                alltips_mask = (1 << ntips) - 1
                class_cached_attr = (frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask)
                type(self)._class_cache[frozenset_leaves] = class_cached_attr
            self._cached_attributes = class_cached_attr
        return self._cached_attributes

    ###############################################################################################

    @property
    def intnodes(self):
        """Lazy evaluation of intnodes"""
        if self._intnodes is None:
            self._intnodes = self.child_dict.keys()
        return self._intnodes

    ###############################################################################################

    @property
    def nodes(self):
        """Lazy evaluation of nodes"""
        if self._nodes is None:
            self._nodes = self.intnodes | self.leaves
        return self._nodes

    ###############################################################################################

    @property
    def frozenset_leaves(self):
        """Return value cached and interned on class.
        Ensures attributes (frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask) are cached and
        interned on class, and copy cached on instance"""

        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes
        return frozenset_leaves

    ###############################################################################################

    @property
    def sorted_leaf_tup(self):
        """Return value cached and interned on class.
        Ensures attributes (frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask) are cached and
        interned on class, and copy cached on instance"""
        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes
        return sorted_leaf_tup

    ###############################################################################################

    @property
    def alltips_mask(self):
        """Return value cached and interned on class.
        Ensures attributes (frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask) are cached and
        interned on class, and copy cached on instance"""
        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes
        return alltips_mask

    ###############################################################################################

    @property
    def leaf2index(self):
        """Return interned value if present, create instance-level cache if not"""
        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes
        return leaf2index

    ###############################################################################################

    @property
    def nodedict(self):
        """Lazy creation of _nodedict when needed"""
        if self._nodedict is None:
            self._nodedict = {node: Nodestruct() for node in self.nodes}
        return self._nodedict

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

        for node in self.leaves:
            remdict[node] = {node}

        for parent in self.sorted_intnodes(deepfirst=False):
            remdict[parent] = set()
            kidstack = list(self.child_dict[parent])
            while kidstack:
                curnode = kidstack.pop()
                if curnode in remdict:
                    remdict[parent].update(remdict[curnode])
                else:
                    kidstack.extend(self.child_dict[curnode])

    ###############################################################################################

    @property
    def remotechildren_indices_dict(self):
        """Dict of {node:frozenset(leaf-indices)}, evaluated lazily when needed.
        Uses remotechildren_dict as starting point"""

        if self._remotechildren_indices_dict is None:
            remidx_dict = self._remotechildren_indices_dict = {}
            getitem = self.leaf2index.__getitem__
            for node, remkids in self.remotechildren_dict.items():
                remidx = frozenset(map(getitem, remkids))
                remidx_dict[node] = remidx
        return self._remotechildren_indices_dict

    ###############################################################################################

    @property
    def remotechildren_mask_dict(self):
        """Dict of {node: remkids_as_bitset}.
        Computed once per tree; fast to query in bipdict() and cladedict()."""

        if self._remotechildren_mask_dict is None:
            frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes

            # remchild[leaf] = leaf, so remchild-masks for leaves = those in leaf2mask
            rembitdict = self._remotechildren_mask_dict = leaf2mask.copy()

            # Mask = union of individual bits in leaves
            # This works since we start iterating from leaves, so child values already computed
            child_dict = self.child_dict
            for parent in self.sorted_intnodes(deepfirst=False):
                mask = 0
                for child in child_dict[parent]:
                    mask |= rembitdict[child]
                rembitdict[parent] = mask

        return self._remotechildren_mask_dict

    ###############################################################################################

    @property
    def rootdist(self):
        """Property (dictionary) giving the distance from each node in tree to the root"""

        if self._rootdist is None:
            self._rootdist = {self.root: 0.0}
            for parent in self.sorted_intnodes(deepfirst=True):
                for child, branch in self.child_dict[parent].items():
                    self._rootdist[child] = (self._rootdist[parent] + branch.length)
        return self._rootdist

    ###############################################################################################

    @property
    def nodedepthdict(self):
        if self._nodedepthdict is None:
            nd = self._nodedepthdict = {}
            rootdist = self.rootdist
            maxdist = 0
            for leaf in self.leaves:
                dist_root_leaf = rootdist[leaf]
                if dist_root_leaf > maxdist:
                    maxdist = dist_root_leaf
            rootdepth = nd[self.root] = maxdist
            for node in self.nodes - {self.root}:
                nd[node] = rootdepth - rootdist[node]
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

    def pathdist_dict(self, rooted=False):
        """Returns dictionary giving path-distance (number of edges) for all pairs of leaves"""

        if rooted and self._pathdist_dict:
            return self._pathdist_dict
        elif not rooted and self._pathdist_dict_unroot:
            return self._pathdist_dict_unroot

        child_dict = self.child_dict
        combinations = itertools.combinations

        # Initialize distdict, counting root traversal as 1 or 2 steps
        if rooted:
            distdict = self._pathdist_dict = {n:{} for n in self.nodes}
            intnodes = self.sorted_intnodes()
        else:
            distdict = self._pathdist_dict_unroot = {n:{} for n in (self.nodes - {self.root})}
            rootkids = self.children(self.root)
            for (child1, child2) in combinations(rootkids, 2):
                distdict[child1][child2] = 1
                distdict[child2][child1] = 1
            intnodes = self.sorted_intnodes()[1:] # Remove root from front of list

        # Build path-distance matrix in optimal order
        # Data structures and algorithm inspired by the Floyd-Warshall algorithm, but modified and
        # faster than O(n^3) since it is on a tree (unique paths)

        # Traverse tree starting from root, breadth-first (sorted_intnodes)
        for parent in intnodes:
            children = child_dict[parent].keys()
            if distdict[parent]:
                prev_contacts = distdict[parent].keys()
                for child in children:
                    for prev_contact in prev_contacts:
                        totlen = distdict[prev_contact][parent] + 1
                        distdict[prev_contact][child] = totlen
                        distdict[child][prev_contact] = totlen
            for (child1, child2) in combinations(children, 2):
                distdict[child1][child2] = 2
                distdict[child2][child1] = 2
            for child in children:
                distdict[parent][child] = 1
                distdict[child][parent] = 1

        return distdict

    ###############################################################################################

    def pathdist_as_ndarray(self, rooted=False):
        """Return flattened version of pathdist_dict (in alphabetical leaf-pair order)"""

        if rooted:
            dictarray = self._pathdist_as_ndarray
        else:
            dictarray = self._pathdist_as_ndarray_unroot
        if dictarray is None:
            distdict = self.pathdist_dict(rooted)
            leafnames = self.sorted_leaf_tup
            namepairs = itertools.combinations(leafnames, 2)
            nleaves = len(leafnames)
            npairs = nleaves * (nleaves - 1) // 2
            distiter = (distdict[n1][n2] for n1, n2 in namepairs)
            dictarray = np.fromiter(distiter, dtype=float, count=npairs)

        return dictarray

    ###############################################################################################

    def copy_treeobject(self, copylengths=True, copyattr=True):
        """Returns copy of Tree object.
        copylengths: copy lengths (otherwise set to 0.0)
        copyattr: copy non-length attributes
        Caches are not copied.
        Similar to effect of copy.deepcopy but customized and much faster"""

        obj = Tree()
        obj.root = self.root
        obj.leaves = self.leaves.copy()
        obj.child_dict = self._copy_child_dict(obj, copylengths, copyattr)
        obj._nodedict = self._copy_nodedict(obj)
        return obj

    ###############################################################################################

    def _copy_child_dict(self, obj, copylengths=True, copyattr=True):
        origdict = self.child_dict
        newdict = {}
        for parent in origdict:
            newdict[parent] = {}
            for child in origdict[parent]:
                origbranch = origdict[parent][child]
                newbranch = origbranch.copy(copylengths, copyattr)
                newdict[parent][child] = newbranch
        return newdict

    ###############################################################################################

    def _copy_nodedict(self, obj):
        if self._nodedict is None:
            return None
        else:
            origdict = self.nodedict
            newdict = {}
            for key,orignode in origdict.items():
                newnode = Nodestruct()
                for attrname, origvalue in orignode.__dict__.items():
                    if isinstance(origvalue, set):
                        newvalue = origvalue.copy()
                    else:
                        newvalue = origvalue  # Assume all other attributes are immutable
                    setattr(newnode, attrname, newvalue)
                newdict[key] = newnode

            return newdict

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
        """Returns sorted intnode list for breadth-first traversal of tree
        deepfirst=False: returns internal nodes in bottom-up (reverse level) order"""

        # "intnodes" is a set, meaning iteration occurs in no defined order.
        # This function returns a list sorted such that deep nodes generally go before
        # shallow nodes (deepfirst=False reverses this)
        cached = self._sorted_intnodes_deep if deepfirst else self._sorted_intnodes_shallow
        if cached is not None:
            return cached

        # Add nodes one tree-level at a time.
        # First root, then children of root, then children of those, etc
        child_dict = self.child_dict
        intnodes = child_dict.keys()

        curlevel = [self.root]
        sorted_nodes = []
        sorted_append = sorted_nodes.append

        while curlevel:
            nextlevel = []
            nextlevel_append = nextlevel.append

            # For each node in current level: add those children that are also internal nodes
            for parent in curlevel:
                sorted_append(parent)
                for child in child_dict[parent]:
                    if child in intnodes:
                        nextlevel_append(child)

            curlevel = nextlevel

        deep = tuple(sorted_nodes)
        shallow = tuple(reversed(sorted_nodes))

        self._sorted_intnodes_deep = deep
        self._sorted_intnodes_shallow = shallow

        return deep if deepfirst else shallow

    ###############################################################################################

    def intnodes_postorder(self):
        """Returns sorted list of intnodes for post-order traversal of tree
        (i.e., a node is only visited after all its children have been visited)"""

        # if self._sorted_intnodes_shallow:
        #     return self._sorted_intnodes_shallow

        out = []
        stack = [(self.root, 0)]  # (node, state) where state=0 means "expand", 1 means "emit"
        while stack:
            node, state = stack.pop()
            if state == 0:
                stack.append((node, 1))
                # push children to expand first
                if node not in self.leaves:
                    for ch in self.children(node):
                        stack.append((ch, 0))
            else:
                if node not in self.leaves:
                    out.append(node)

        # self._sorted_intnodes_shallow = out
        return out

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
        """Returns dict_keys object (similar to set) containing parent's immediate descendants"""

        # Python note: does not seem to benefit from lru_caching, and leads to multiple problems

        try:
            return self.child_dict[parent].keys()
        except KeyError as err:
            raise TreeError(f"Node {parent} is not an internal node") from err

    ###############################################################################################

    def child_edges(self, parent):
        """Returns child_edges dict {kid1:branchstruct, kid2:branchstruct [...]} for parent"""

        return self.child_dict[parent]

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

    def is_parent_child_pair(self, pnode, cnode):
        """Returns True if cnode is child of pnode. False otherwise (also if order reversed)"""

        # Python note: hack. There must be a better solution (or naming, or call sequence) to this.
        # Or maybe I should never have to ask...

        try:
            if pnode == self.parent(cnode):
                return True
            else:
                return False
        except TreeError:
            return False


    ###############################################################################################

    def branch_set(self):
        """Returns set of (parent, child) tuples for all branches in tree"""

        branch_set = set()
        for p in self.intnodes:
            for c in self.children(p):
                branch_set.add((p,c))
        return branch_set

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
            self = self.copy_treeobject(copylengths=False, copyattr=False)
            other = other.copy_treeobject(copylengths=False, copyattr=False)

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

        # special case: mrca for leaf is leaf itself
        if len(leaves) == 1:
            return next(iter(leaves))

        # python note: is if test always faster than set conversion?
        if isinstance(leaves, set):
            leafset = leaves
        else:
            leafset = set(leaves)
        remotechildren_dict = self.remotechildren_dict
        parent_dict = self.parent_dict

        # pick random starting node among leafset, and find its parent node
        random_leaf = next(iter(leafset))
        parent = parent_dict[random_leaf]

        # Walk down the tree from the initially picked node, until remkids include all of "leafset"
        while not leafset <= remotechildren_dict[parent]:
            parent = parent_dict[parent]

        return parent

    ###############################################################################################

    def find_mrca_mask(self, query_mask, start_leaf):
        """MRCA using bitmasks.

        Relies on caller having precomputed query_mask using same sorted_leaf_tup as self

        query_mask: set of leaves represented as a bitmask, i.e., an int with bits set on
                    relevant positions (position = index in sorted leaf tuple)

        start_leaf: which leaf in query to start climbing towards root from (random)"""

        remmask = self.remotechildren_mask_dict
        parent_dict = self.parent_dict
        parent = parent_dict[start_leaf]
        while (remmask[parent] & query_mask) != query_mask:  # Is query_mask subset?
            parent = parent_dict[parent]
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

        leaves1 = set(leafset)
        if leaves1 == self.leaves:
            return self.root
        else:
            leaves2 = self.leaves - leaves1
            base1, base2 = self.find_bipart_nodes(leaves1, leaves2, check_arguments=True)
            return base1

    ###############################################################################################

    def find_bipart_nodes(self, leaves1, leaves2, check_arguments=True):
        """Find two internal nodes corresponding to given bipartition of leaves.

        Input: leaves1, leaves2: leaves on either side of internal branch (bipartition)
        Out: (basenode1, basenode2) nodes on either end of the branch corresponding to the bipartition
             basenode1 is at the base of leaves1 leaves
             basenode2 is at the base of leaves2 leaves

        check_arguments: check that leaves1, leaves2 are valid bipartition of leaves

        If tree is rooted at bifurcation, and bipartition corresponds to two halves of tree,
        then the two kids of root will be the basenodes (both basal branches are considered
        to be part of the same branch)"""

        if check_arguments:
            bip1 = set(leaves1)
            bip2 = set(leaves2)
            all_leaves = bip1 | bip2
            if all_leaves < self.leaves:
                missing_leaves = self.leaves - all_leaves
                raise TreeError(f"Not a bipartition: some leaves are not included in either bipartition1 or 2: {missing_leaves}")
            if self.leaves < all_leaves:
                superfluous = all_leaves - self.leaves
                raise TreeError(f"Not a bipartition: these leaves are not on tree: {superfluous}")
            if not bip1:
                raise TreeError("Not a bipartition: leaves1 contains no leaves")
            if not bip2:
                raise TreeError("Not a bipartition: leaves2 contains no leaves")

        mrca1 = self.find_mrca(leaves1)
        mrca2 = self.find_mrca(leaves2)
        if mrca1 != self.root and mrca2 == self.root:
            return mrca1, self.parent(mrca1)
        elif mrca1 == self.root and mrca2 != self.root:
            return self.parent(mrca2), mrca2
        else:
            return mrca1, mrca2   # Both are kids of root or mrca1==mrca2 (bip anchored in multifurcation)

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
        return cumdist2 + node1_ancdist[child2]

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

        return self.nodedepthdict[node]

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

        # Find path from node1 to root
        root = self.root
        parent_dict = self.parent_dict
        path1 = [node1]
        curnode = node1
        while curnode != root:
            curnode = parent_dict[curnode]
            path1.append(curnode)
        path1set = set(path1)

        # Find path from node2 to root (or to first node that is also on node1's path)
        path2 = [node2]
        curnode = node2
        while curnode not in path1set:
            curnode = parent_dict[curnode]
            path2.append(curnode)

        # merge paths
        intersect = path1.index(curnode)
        path2.reverse()
        fullpath =  path1[:intersect] + path2

        return fullpath

    ###############################################################################################

    def find_path_intersection(self, node1, node2):
        """Returns intersection between paths from node1 and node2 to root (two-leaf mrca)."""

        # Python note: overlap with nodepath - integrate somehow perhaos?

        # Find nodes on path from node1 to root
        root = self.root
        parent_dict = self.parent_dict
        path1set = set([node1])
        curnode = node1
        while curnode != root:
            curnode = parent_dict[curnode]
            path1set.add(curnode)

        # Find path from node2 to root (or to first node that is also on node1's path)
        curnode = node2
        while curnode not in path1set:
            curnode = parent_dict[curnode]

        return curnode

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

        # Python note: check for zero branch lengths somehow
        # using .length() is a bit expensive. Some sort of try except approach?

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

        # Construct dictionary mapping old to new names (new names = shuffled version of old)
        oldnames = list(self.leaves)
        newnames = list(self.leaves)
        random.shuffle(newnames)
        old2new = dict(zip(oldnames, newnames))

        # Loop through child_dict and replace leafnames where they occur
        for parent in self.child_dict:
            tmp_dict = {}
            for child, value in self.child_dict[parent].items():
                if child in self.leaves:
                    tmp_dict[old2new[child]] = value
                else:
                    tmp_dict[child] = value
            self.child_dict[parent] = tmp_dict

        # Delete all caches which will now be rendered obsolete
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

    def newick(self, printdist=True, printlabels=True, labelfield="label", precision=6,
               transdict=None, node_attributes=None, branch_attributes=None):
        """Returns Newick format tree string representation of tree object, with optional metacomments"""

        def create_metacomment(struct, attributes):
            """Helper function to create metacomment strings based on attributes of a structure"""
            tmplist = []
            for attrname in attributes:
                try:
                    value = getattr(struct, attrname)
                except AttributeError:
                    msg = (f"'{struct.__class__.__name__}' object has no attribute '{attrname}'.\n"
                          f"Available attributes: {list(vars(struct).keys())}")
                    raise TreeError(msg)
                if isinstance(value, float):
                    tmplist.append(f"{attrname}={value:.{precision}g}")
                else:
                    tmplist.append(f"{attrname}={value}")
            return f"[&{','.join(tmplist)}]"

        def append_children(parentnode):
            """Recursive function that has main responsibility for building Newick tree string"""

            for child in self.children(parentnode):
                branchstruct = self.child_dict[parentnode][child]
                dist = branchstruct.length
                if child in self.leaves:
                    treelist.append(transdict[child] if transdict else child)
                else:
                    treelist.append("(")
                    append_children(child)
                    treelist.append(")")
                    if printlabels:
                        label = getattr(branchstruct, labelfield, "")
                        if isinstance(label, float):
                            treelist.append(f"{label:.{precision}g}")
                        else:
                            treelist.append(f"{label}")
                if node_attributes:
                    metacomment_node = create_metacomment(self.nodedict[child], node_attributes)
                    treelist.append(metacomment_node)
                if printdist:
                    treelist.append(f":{dist:.{precision}g}")
                if branch_attributes:
                    metacomment_branch = create_metacomment(branchstruct, branch_attributes)
                    treelist.append(metacomment_branch)

                treelist.append(",")
            del treelist[-1]  # Remove last comma when no more siblings

        # EXECUTION STARTS HERE!
        if node_attributes and (self._nodedict is None):
            raise TreeError(f"Tree has no nodedict. Can not print node-related attributes: {node_attributes}")
        root = self.root
        treelist = ["("]
        append_children(root)
        if node_attributes:
            metacomment_node = create_metacomment(self.nodedict[root], node_attributes)
            treelist.append(f"){metacomment_node};")
        else:
            treelist.append(");")
        treestring = "".join(treelist)
        return treestring

    ###############################################################################################

    def nexus(self, printdist=True, printlabels=True, labelfield="label", precision=6,
              translateblock=False,  node_attributes=None, branch_attributes=None,
              colorlist=None, colorfg="#0000FF", colorbg="#000000"):
        """Returns nexus format tree as a string"""

        # Construct header
        if colorlist:
            header = "#NEXUS\nbegin taxa\n\tdimensions ntax={};\n".format(len(self.leaflist()))
            header += "\ttaxlabels\n\t"
            stringlist = [header]
            for leaf in self.leaflist():
                stringlist.append("\t{}".format(leaf))
                col = colorfg if leaf in colorlist else colorbg  #
                stringlist.append(f"[&!color={col}]\n")
            stringlist.append(";\nend;\n")
        else:
            stringlist = ["#NEXUS\n\nbegin trees;\n"]

        # If translateblock is requested: add translateblock to stringlist
        transdict = None
        if translateblock:
            transdict = self.transdict()
            stringlist.append(self.translateblock(transdict))

        # Add newick tree string with optional meta-comments for figtree format
        stringlist.append("\ttree nexus_tree = ")
        stringlist.append(self.newick(printdist=printdist, printlabels=printlabels, labelfield=labelfield,
                                      precision=precision, transdict=transdict,
                                      node_attributes=node_attributes, branch_attributes=branch_attributes))

        # Add footer
        stringlist.append("\nend;\n")

        return "".join(stringlist)

    ###############################################################################################

    def iter_bipinfo(self, keep_remchild_dict=False):
        """
        Stream bipartitions for all branches.
        Yields: (bipartition, branchstruct)
        Avoids building a bipdict.
        """
        remmask = self.remotechildren_mask_dict
        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes
        from_halfmask = Bipartition.from_halfmask
        child_dict = self.child_dict

        object_cache = Bipartition._class_cache.get(frozenset_leaves)
        if object_cache is None:
            object_cache = {}
            Bipartition._class_cache[frozenset_leaves] = object_cache

        root = self.root
        rootkids = self.children(root)

        # If exactly two root kids, we will skip yielding the two root edges and yield the merged one later
        do_root_merge = (len(rootkids) == 2)
        if do_root_merge:
            kid1, kid2 = rootkids

        for parent, child_edges in child_dict.items():
            # skip the two root edges (root->kid1 and root->kid2) if we will merge them
            if parent == root and do_root_merge:
                continue

            for child, branch in child_edges.items():
                halfmask = remmask[child]
                bip = from_halfmask(halfmask, alltips_mask, object_cache, sorted_leaf_tup)
                yield bip, branch

        if do_root_merge:
            rootbip = self.rootbip()
            branch1 = child_dict[root][kid1]
            branch2 = child_dict[root][kid2]
            branch_merged = branch1.merge(branch2, check_compat=True)
            yield rootbip, branch_merged

        if not keep_remchild_dict:
            self._remotechildren_mask_dict = None

    ###############################################################################################

    def bipdict(self, keep_remchild_dict = False):
        """Deprecated: use generator/streaming version iter_bipinfo when possible.
        Returns tree in the form of a "bipartition dictionary"
        """

        return dict(self.iter_bipinfo(keep_remchild_dict=keep_remchild_dict))

    ###############################################################################################

    def iter_cladeinfo(self, node2clade=None, keep_remchild_dict=False):
        """
        Stream clade info for each node.

        Yields:
            (node, clade, depth, nleaves)

        Side effects:
            If node2clade is a dict, fills node2clade[node] = clade.

       This avoids building cladedict and avoids allocating Nodestruct per node per tree
        """

        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes
        from_mask = Clade.from_mask

        object_cache = Clade._class_cache.get(frozenset_leaves)
        if object_cache is None:
            object_cache = {}
            Clade._class_cache[frozenset_leaves] = object_cache

        # Local binds for speed
        remkid_items = self.remotechildren_mask_dict.items
        nodedepth = self.nodedepth

        if node2clade is None:
            for node, remkids_mask in remkid_items():
                clade = from_mask(remkids_mask, object_cache, sorted_leaf_tup)
                yield node, clade, nodedepth(node), remkids_mask.bit_count()
        else:
            # fill node2clade
            for node, remkids_mask in remkid_items():
                clade = from_mask(remkids_mask, object_cache, sorted_leaf_tup)
                node2clade[node] = clade
                yield node, clade, nodedepth(node), remkids_mask.bit_count()

        if not keep_remchild_dict:
            self._remotechildren_mask_dict = None

    ###############################################################################################

    def cladedict(self, keep_remchild_dict=False, track_subcladepairs=False):
        """
        Deprecated: use generator version iter_cladeinfo when possible.
        This function only to maintain backwards-compatible materialization of dict

        Returns tree in the form of a "clade dictionary": {clade: nodestruct}
        Nodestructs get the .depth and ntips attributes set during construction.
        If requested: Nodestructs also sets .subcladepairs attribute (set)
        """
        cladedict = {}
        node2clade = {} if track_subcladepairs else None
        it = self.iter_cladeinfo(node2clade=node2clade, keep_remchild_dict=keep_remchild_dict)

        Nodestruct_ = Nodestruct  # local bind
        for node, clade, depth, nleaves in it:
            cladedict[clade] = Nodestruct_(depth, nleaves)

        if track_subcladepairs:
            # Second pass: attach observed child clade pairs to the per-tree nodestructs
            child_dict = self.child_dict
            intnodes = self.intnodes

            for node in intnodes:
                try:
                    kid1, kid2 = child_dict[node].keys()
                except ValueError:
                    msg = (
                        f"Error while tracking subclade_pairs - clade has more than 2 children: "
                        f"parent-node: {node}    children: {self.children(node)}"
                        f"Are you sure input trees are rooted, e.g. using a clock-model?"
                    )
                    raise TreeError(msg)

                parent_clade = node2clade[node]
                ns = cladedict[parent_clade]
                if ns.nleaves > 2:
                    ns.add_subcladepair(node2clade[kid1], node2clade[kid2])

        return cladedict

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

    def _rootbip_core(self, trackrootblen=False):
        """Helper method for rootbip and rootbip_with_frac.
        
        Returns 
            rootbip, root-kid1, root-kid2, halfmask
        
        where halfmask is the canonical halfmask for root bipartition
        """
        
        root = self.root
        rootkids = self.child_dict[root]
        if len(rootkids) != 2:
            raise TreeError(
                "Tree rooted at multifurcation - not possible to assign a unique root-bipartition")

        kid1, kid2 = rootkids.keys()
        halfmask = self.remotechildren_mask_dict[kid1]
        
        frozenset_leaves, sorted_leaf_tup, leaf2index, leaf2mask, ntips, alltips_mask = self.cached_attributes
        object_cache = Bipartition._class_cache.get(frozenset_leaves)
        if object_cache is None:
            object_cache = {}
            Bipartition._class_cache[frozenset_leaves] = object_cache        
        rootbip = Bipartition.from_halfmask(halfmask, alltips_mask, object_cache, sorted_leaf_tup)

        return rootbip, kid1, kid2, halfmask
            
    ###############################################################################################

    def rootbip(self):
        """For a tree rooted at a bifurcation: returns Bipartition object corresponding to 
        bipartition present at root"""
        
        rootbip, kid1, kid2, halfmask = self._rootbip_core()
        return rootbip

    ###############################################################################################

    def rootbip_with_frac(self):
        """For a tree rooted at a bifurcation: returns Bipartition object and blen info
                (Bipartition, frac_to_canon)
        frac_to_canon = blen_canonical / (blen1 + blen2)
        where "canonical" means the side whose mask equals rootbip._mask
        """
        
        rootbip, kid1, kid2, halfmask = self._rootbip_core()

        root_child_edges = self.child_dict[self.root]
        blen1 = root_child_edges[kid1].length
        blen2 = root_child_edges[kid2].length

        denom = blen1 + blen2
        frac1 = 0.5 if denom <= 0.0 else blen1 / denom

        # canonical side = rootbip._mask
        frac_to_canon = frac1 if (halfmask == rootbip._mask) else (1.0 - frac1)

        return rootbip, frac_to_canon

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
            # Note: not safe to use .children() method while changing tree (cache will break)
            numkids = len(self.child_dict[node])
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
                intnode2 = self.insert_node(intnode1, subset1, Branchstruct())
                if subset1_size > 2:
                    unresolved_nodes.append(intnode2)

            if subset2_size > 1:
                intnode2 = self.insert_node(intnode1, subset2, Branchstruct())
                if subset2_size > 2:
                    unresolved_nodes.append(intnode2)

    ###############################################################################################

    def set_branch_attribute(self, node1, node2, attrname, attrvalue):
        """Set the value of any branch attribute.
        attrname: Name of attribute (e.g., "length")
        attrvalue: Value of attribute (e.g. 0.153)"""

        setattr(self.child_dict[node1][node2], attrname, attrvalue)
        self.clear_caches()

    ###############################################################################################

    def get_branch_attribute(self, node1, node2, attrname, default=""):
        """Get the value of any branch attribute.
        attrname: Name of attribute (e.g., "length")"""

        return getattr(self.child_dict[node1][node2], attrname, default)

    ###############################################################################################

    def set_node_attribute(self, node, attrname, attrvalue):
        """Set attribute for the specified node.
        attrname: Name of attribute (e.g., "depth")
        attrvalue: Value of attribute (e.g. 0.153)"""

        setattr(self.nodedict[node], attrname, attrvalue)

    ###############################################################################################

    def get_node_attribute(self, node, attrname, default=""):
        """Get specified attribute for the specified node"""

        return getattr(self.nodedict[node], attrname, default)

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
        self.clear_length_caches()

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

        return getattr(self.child_dict[parent][child], "label", "")

    ###############################################################################################

    def set_nodeid_labels(self):
        """Sets labels to be the same as the child node ID
        Allows use of e.g. Figtree to show nodeIDs as nodelabels"""

        for parent in self.intnodes:
            for kid in self.children(parent):
                self.setlabel(parent, kid, str(kid))

    ###############################################################################################

    def set_blens_from_depths(self):
        """Set all branch lengths based on node depths: blen = depth_parent - depth_child"""

        # Ensure presence of nodedict and that depth attribute exists for all nodes
        if self._nodedict is None:
            raise TreeError("Tree does not have nodedict, so depth information not available")
        missing = [n for n in self.nodes if not hasattr(self.nodedict[n], "depth")]
        if missing:
            msg = (f"No information about node depths on this tree "
                   f"(missing 'depth' attribute for {len(missing)} nodes; e.g. {missing[:5]})")
            raise TreeError(msg)

        # Set lengths from the stored depths
        for parent in self.sorted_intnodes(deepfirst=True):
            p_depth = self.nodedict[parent].depth
            for child in self.children(parent):
                c_depth = self.nodedict[child].depth
                self.child_dict[parent][child].length = (p_depth - c_depth)
        self.clear_length_caches()

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
            other.root = basenode
            curlevel = [basenode]
            while curlevel:
                nextlevel = []
                for parent in curlevel:
                    other.child_dict[parent] = {}
                    kids = self.children(parent)
                    for kid in kids:
                        other.child_dict[parent][kid] = self.child_dict[parent][kid].copy()
                    intnode_kids = kids & self.intnodes
                    nextlevel.extend(intnode_kids)
                    leaf_kids = kids & self.leaves
                    other.leaves.update(leaf_kids)
                curlevel = nextlevel

        # If self.nodedict exists: copy relevant parts from self to other
        if self._nodedict is not None:
            for node in other.nodes:
                other.nodedict[node] = self.nodedict[node]

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

        # If node2 is not root of tree2 (or child of root2) then re-root on branch below node2.
        # After this, grafting can happen below root2
        elif (node2 != other.root) and (node2 not in other.children(other.root)):
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

        # Update set of leaves
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

        # Account for fact that some or all intnodes may be strings (e.g., transmission trees)
        if all(isinstance(intnode, str) for intnode in self.intnodes):
            newnode = 0
        else:
            newnode = max(x for x in self.intnodes if isinstance(x, int)) + 1
        tree[newnode] = {}

        # Add new internal node as child of "parent"
        tree[parent][newnode] = branchstruct

        # Move childnodes from previous parent to new node
        for child in childnodes:
            tree[newnode][child] = tree[parent][child]
            del tree[parent][child]

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
            branchstruct1 = branchstruct.copy()
            branchstruct1.length = branchstruct.length / 2
            branchstruct2 = branchstruct.copy()
            branchstruct2.length = branchstruct.length / 2
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
            newbranch = self.child_dict[child][grandchild].copy()
            newbranch.length += addlen
            self.child_dict[parent][grandchild] = newbranch

        # Delete "child" node and link from parent to child.
        del self.child_dict[child]
        del self.child_dict[parent][child]

        self.clear_caches()

    ###############################################################################################

    def remove_leaves(self, leaflist):
        """Removes leaves in list from tree, cleans up so remaining tree structure is sane"""

        for leaf in leaflist:
            self.remove_leaf(leaf)

    ###############################################################################################

    def remove_leaf(self, leaf):
        """Removes named leaf from tree, cleans up so remaining tree structure is sane.
        Internal nodes may be removed, but will not change their IDs
        Note: assumes no unary nodes (leaf has at least one sibling)"""

        parent = self.parent(leaf)
        childset = self.children(parent)
        root = self.root

        if self._nodedict is not None:
            orignodes = self.nodes.copy()

        # If leaf is part of bifurcation AND is directly attached to root, then
        # the "other child" of the root must become the new root
        if (len(childset) == 2) and (leaf in self.children(root)):
            [child2] = childset - {leaf}                    # Remaining item is other child
            del self.child_dict[root]                       # Remove entry for old root
            self.root = child2                              # child2 is new root

        # If leaf is part of bifurcation but NOT attached directly to root, then parent
        # must also be removed from tree, and the remaining child needs to be grafted
        # onto grandparent with proper cumulated distance
        # Keep branch label (if any) of the branch from grandparent to parent
        elif len(childset) == 2:
            [child2] = childset - {leaf}                    # Remaining item is other child
            child2dist = self.child_dict[parent][child2].length   # Remember dist to other child
            grandparent = self.parent(parent)

            # Add remaining child to grandparent
            self.child_dict[grandparent][child2] = self.child_dict[grandparent][parent]
            self.child_dict[grandparent][child2].length += child2dist   # Cumulated distance
            del self.child_dict[parent]                      # Delete parent and leaf
            del self.child_dict[grandparent][parent]         # Also remove pointer from gp to p
            del self._parent_dict[leaf]                      # Remove unused entries in parent_dict
            del self._parent_dict[parent]

        # If leaf is part of multifurcation, then no special cleanup needed
        else:
            del self.child_dict[parent][leaf]
            del self._parent_dict[leaf]

        # Remove leaf entry from global leaflist. Update intnodeslist
        self.leaves.remove(leaf)

        # Clean up nodedict if present
        if self._nodedict is not None:
            remnodes = orignodes - self.nodes
            for node in remnodes:
                del self.nodedict[node]

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

    def nodepath_to_branchset(self, nodepath):
        """Input: nodepath (list of (node1, node2) tuples from .nodepath method)
        Output: set of branches (set of (parent, child) tuples)"""

        branchset = set()
        for i in range( len(nodepath) - 1 ):
            n1, n2 = nodepath[i], nodepath[i+1]
            if self.is_parent_child_pair(n1, n2):
                branchset.add((n1,n2))
            else:
                branchset.add((n2,n1))
        return branchset

    ###############################################################################################

    def _init_with_keeplist(self, keeplist):
        """Helper function for prune_maxlen. Initilalizes used_branches, keep_leaves when
        keeplist provided"""
        used_branches = set()

        # Special case: only 1 leaf in keeplist
        # find most distant second leaf and add to keeplist. Then algorithm works as usual
        if len(keeplist) == 1:
            leaf1 = keeplist[0]
            leaf2 = self.find_most_distant(leaf1, self.leaves - {leaf1})[0]
            keeplist.append(leaf2)

        # Find mrca and add all branches on spread out subtree to used_branches
        keeplist_mrca = self.find_mrca(keeplist)
        for leaf in keeplist:
            path = self.nodepath(leaf, keeplist_mrca)
            path_branches = self.nodepath_to_branchset(path)
            used_branches |= path_branches
        keep_leaves = set(keeplist)
        return used_branches, keep_leaves

    ###############################################################################################

    def _init_with_longestbranch(self):
        """Helper function for prune_maxlen: Initializes used_branches, keep_leaves based on
        finding longest path. Note: longest path is the optimal pruned tree for 2 leaves"""

        maxdist,node1,node2 = self.diameter(return_leaves=True)
        maxpath = self.nodepath(node1, node2)
        used_branches = self.nodepath_to_branchset(maxpath)
        keep_leaves = set([node1, node2])
        return used_branches, keep_leaves

    ###############################################################################################

    def _init_possible(self, used_branches, keep_leaves_mrca):
        """Helper function for prune_maxlen: initialize possible_basal_branches based on
        used_branches and keep_leaves.
        Possible basal branches includes branches sprouting off used_branches.
        If root is not already in used_branches, then possible_basal_branches also includes
        branches on path to root and sprouting from that path"""

        possible_basal_branches = set()
        for (p,c1) in used_branches:
            for c2 in (self.children(p) - {c1}):
                if (p,c2) not in used_branches:    # Python note: faster to just set subtract at end?
                    possible_basal_branches.add((p,c2))

        # If root not in used_branches: add path to root, and branches sprouting from this path
        if keep_leaves_mrca != self.root:
            rootpath = self.nodepath(self.root, keep_leaves_mrca)
            rootpath_branches = self.nodepath_to_branchset(rootpath)
            possible_basal_branches.update(rootpath_branches)
            for (p,c1) in rootpath_branches:
                for c2 in (self.children(p) - {c1}):
                    if (p,c2) not in used_branches:
                        possible_basal_branches.add((p,c2))

        return possible_basal_branches

    ###############################################################################################

    def _find_side_branches(self, branchset):
        """Helper function for prune_maxlen. Finds all side branches to provided set of branches"""
        side_branches = set()
        for p,c1 in branchset:
            for c2 in (self.children(p) - {c1}):
                side_branches.add((p,c2))
        return side_branches - branchset

    ###############################################################################################

    def _find_possible_down(self, keep_leaves_mrca):
        """Helper function for prune_maxlen. Initializes possible_up_basal and possible_down_basal"""

        root_path = self.nodepath(self.root, keep_leaves_mrca)
        root_path_branches = self.nodepath_to_branchset(root_path)
        possible_down_basal = self._find_side_branches(root_path_branches)

        return possible_down_basal

    ###############################################################################################

    def _check_prune_maxlen_arguments(self, nkeep, keeplist):
        """Helper function for prune_maxlen. Checks sanity of arguments, raises error if issue"""

        if nkeep > len(self.leaves):
            raise TreeError(f"nkeep > number of leaves: {nkeep} > {len(self.leaves)}")
        if len(keeplist) > len(self.leaves):
            raise TreeError(f"len(keeplist) > number of leaves: {len(keeplist)} > {len(self.leaves)}")
        if len(keeplist) > nkeep:
            raise TreeError( f"len(keeplist) > nkeep: {len(keeplist)} > {nkeep}")

    ###############################################################################################

    def prune_maxlen(self, nkeep, keeplist=[], return_leaves=False):
        """Prune tree so remaining nkeep leaves spread out maximal percentage of branch length
        keeplist: optional list of leaves that must be included.
        return_leaves: return selected leaves without pruning tree
        Note: Best solution including keeplist may be less good than optimal solution"""

        self._check_prune_maxlen_arguments(nkeep, keeplist)

        # Initialize
        if keeplist:
            used_branches, keepleaves = self._init_with_keeplist(keeplist)
        else:
            used_branches, keepleaves = self._init_with_longestbranch()
        keepleaves_mrca = self.find_mrca(keepleaves)
        possible_up_basal = self._find_side_branches(used_branches)
        possible_down_basal = self._find_possible_down(keepleaves_mrca)

        # Until we have added nkeep leaves to path:
        # find longest newpath from existing path to leaf, add to path, update vars
        while len(keepleaves) < nkeep:
            keepleaves_mrca_changed = False
            maxdist = 0.0
            for parent,child in possible_up_basal:
                for leaf in self.remote_children(child):
                    dist = self.nodedist(parent, leaf)
                    if dist > maxdist:
                        keep_parent, keep_leaf, maxdist = parent, leaf, dist
            for parent,child in possible_down_basal:
                downdist = self.nodedist(parent, keepleaves_mrca)
                for leaf in self.remote_children(child):
                    dist = self.nodedist(parent, leaf) + downdist
                    if dist > maxdist:
                        keep_parent, keep_leaf, maxdist = keepleaves_mrca, leaf, dist
                        keepleaves_mrca_changed = True
            keepleaves.add( keep_leaf )

            # Update relevant variables
            if keepleaves_mrca_changed:
                keepleaves_mrca = self.find_mrca(keepleaves)
                possible_down_basal = self._find_possible_down(keepleaves_mrca)
            newpath = self.nodepath(keep_parent, keep_leaf)
            newpath_branches = self.nodepath_to_branchset(newpath)
            used_branches |= newpath_branches
            newpath_sidebranches = self._find_side_branches(newpath_branches) - used_branches
            possible_up_basal |= newpath_sidebranches
            possible_up_basal -= used_branches

        self.clear_caches()

        # If requested: return selected leaves without pruning tree
        # Otherwise: prune tree so only leaves in keepset are retained
        if return_leaves:
            return keepleaves
        else:
            discardleaves = self.leaves - keepleaves
            self.remove_leaves(discardleaves)

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
        self.child_dict[parent][newname] = self.child_dict[parent][oldname]
        del self.child_dict[parent][oldname]

        # Update self.leaves
        self.leaves.add(newname)
        self.leaves.remove(oldname)

        # Update self.parent_dict
        self.parent_dict[newname] = self.parent_dict[oldname]
        del self.parent_dict[oldname]

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

    def treedist_pathdiff(self, other, rooted=False):
        """Compute path difference tree-distance between self and other:
        Euclidean distance between nodepath-dist matrices considered as vectors.
        Measure described in M.A. Steel, D. Penny, Syst. Biol. 42 (1993) 126–141
        rooted=True: count traversal of root as 2 steps (not 1)
        """

        self_distvec = self.pathdist_as_ndarray(rooted)
        other_distvec = other.pathdist_as_ndarray(rooted)

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
            branch1 = self.child_dict[root][kid1]
            branch2 = self.child_dict[root][kid2]
            branch_merged = branch1.merge(branch2, check_compat=True)   # Sums up branch lengths, merges attributes
            if kid1 in self.intnodes:
                self.child_dict[kid1][kid2] = branch_merged
                self.root = kid1
            elif kid2 in self.intnodes:
                self.child_dict[kid2][kid1] = branch_merged
                self.root = kid2
            else:
                raise TreeError("Cannot deroot tree: only leaves are present")

            # remove previous root
            del self.child_dict[root]

            # clear caches, which are now unreliable
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
                msg = "Need to specify node2 to reroot() method when not rooting at polytomy"
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
            newbranch = self.child_dict[parent][child].copy()
            newbranch.length = parent_to_root_dist
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
        leaves = tuple(self.leaves)
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
        minchild = next(iter(self.children(minparent)))
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
                    rootkids = rootkids - {minchild}
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
        if isinstance(outgroup, str):
            outgroup = set([outgroup])
        else:
            outgroup = set(outgroup)
        ingroup = self.leaves - outgroup
        inbase, outbase = self.find_bipart_nodes(ingroup, outgroup)

        # If outgroup should form basal polytomy with ingroup: root on outbase
        if polytomy:
            self.reroot(outbase, inbase, polytomy=True)

        # If tree does not have branch lengths: skip computation of lengths
        elif self.length() == 0.0:
            self.reroot(outbase, inbase)

        # Else: place root in the middle between inbase and outbase
        else:
            half_dist = self.nodedist(outbase, inbase) / 2
            self.reroot(outbase, inbase, node1dist=half_dist)

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
        possible_regraft_nodes = treecopy.nodes - {treecopy.root}
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

        # Select random regraft_node or check the one provided
        possible_regraft_nodes = self.possible_spr_regraft_nodes(prune_node)
        if regraft_node is None:
            regraft_node = random.choice(list(possible_regraft_nodes))
        elif regraft_node not in possible_regraft_nodes:
            msg = f"Specified regraft_node {regraft_node} is not compatible with prune_node"
            raise TreeError(msg)

        # Pruning: Remove subtree
        isleaf = prune_node in self.leaves      # Has to be set before pruning!
        subtree = self.prune_subtree(prune_node)

        # Regraft: Add subtree back onto remaining tree
        # Special treatment when pruning single leaf (to avoid superfluous internal node)
        self.graft(subtree, regraft_node, graft_with_other_root=isleaf)

    ###############################################################################################

    def prune_subtree(self, basenode):
        """Prune subtree rooted at basenode from self. Returns pruned subtree"""

        subtree = self.subtree(basenode)
        for leaf in self.remote_children(basenode):
            self.remove_leaf(leaf)
        return subtree

    ###############################################################################################

    def parsimony_possible_states(self):
        """Performs parsimony analysis on tree object, and returns info about possible states
        at internal nodes.

        Assumes tree object has .nodedict attribute: {nodeid:Nodestruct}, and that each
        Nodestruct has a .state attribute, which is either empty (unlabelled nodes)
        or a string specifying the state (e.g., "mink" or "human" or "A").

        Returns: copy of tree object where nodestructs now also have attributes
        .primary_set, .secondary_set, and .optimal_set
        .optimal_set gives the full set of states that a node can have on the MP tree.
        Note: not all combinations of states will yield maximum parsimony score.

        Based on Hartigan's 1973 parsimony algorithm.
        This version allows multifurcations, and internal nodes can have observed state.

        Note: so far only implemented for one-character states (in the cladistic sense).
        """

        # Note: i am using algorithm from:
        # Live phylogeny with polytomies: Finding the most compact parsimonious trees
        # D. Papamichaila et al., Computational Biology and Chemistry 69 (2017) 171–177
        # Upper and lower sets are here referred to as primary and secondary set

        if self._nodedict is None:
            msg = "Tree object has no .nodedict attribute. Can't perform parsimony analysis"
            raise TreeError(msg)

        treecopy = self.copy_treeobject()
        treecopy = self._hartigan_tip2root_pass(treecopy)
        treecopy = self._hartigan_root2tip_pass(treecopy)

        return treecopy

    ###############################################################################################

    def parsimony_assign_fits(self, fitpref=None):
        """Performs parsimony analysis, and returns copy of tree object where one optimal fit
        (i.e., one giving the maximum parsimony number of changes) has been assigned to each node.
        fitpref = None: choose random fits at unlabelled nodes, where possible
        fitpref = <concrete state-string>: Set ambiguous state at unlabelled nodes to this value
                 when possible (i.e., when compatible with MP score)
        """

        tree = self.parsimony_possible_states()
        ndict = tree.nodedict

        # Set fit and wasambig for labelled nodes
        for node in tree.nodes:
            if ndict[node].state:
                ndict[node].fit = ndict[node].state
                ndict[node].wasambig = False

        # Set fit and wasambig for root
        self._set_fit_and_ambiguity_for_node(tree, tree.root, None, fitpref)

        # Set fit and wasambig for internal, unlabelled nodes
        for parent in tree.sorted_intnodes(deepfirst=True):
            for child in tree.children(parent):
                if not ndict[child].state:
                    self._set_fit_and_ambiguity_for_node(tree, child, parent, fitpref)

        return tree

    ###############################################################################################

    def parsimony_count_changes(self, fitpref=None):
        """Performs parsimony analysis on tree object, and returns parsimony score and
        dictionary giving number of changes in either direction:
            key = tuple (from,to), value = count in that direction
            {(state1, state2):count1_to_2,
             (state2, state1):count2_to_1}

        fitpref = None: choose random fits at unlabelled nodes, where possible
        fitpref = <concrete state-string>: Set ambiguous state at unlabelled nodes to this value
                 when possible (i.e., when compatible with MP score)
        """

        tree = self.parsimony_assign_fits(fitpref)
        pscore, countdict = self._count_changes_on_fitted_tree(tree)
        return pscore, countdict

    ###############################################################################################

    def _count_changes_on_fitted_tree(self, tree):
        """Helper function for parsimony analyses:
        Input: tree where nodedict already has assigned fits for each node
        Output: tuple of parsimony-score and countdict. {(from_state, to_state):count, ...}
        key = tuple of (from_state,to_stat), value = count in that direction"""

        ndict = tree.nodedict
        countdict = defaultdict(int)
        pscore = 0
        for parent in tree.sorted_intnodes(deepfirst=True):
            pfit = ndict[parent].fit
            for child in tree.children(parent):
                cfit = ndict[child].fit
                if pfit != cfit:
                    countdict[(pfit,cfit)] += 1
                    pscore += 1

        return pscore, countdict

    ###############################################################################################

    def _hartigan_tip2root_pass(self, tree):
        """Perform tip-to-root pass of Hartigan parsimony algorithm.
        Input: tree object with .nodedict having .state info for each node (some empty)
        Output: same tree object with added attribute .primary_set for all nodes
        """

        # Set primary set of all labelled nodes (both internal and leafs) = state-value
        # Set secondary set of labelled nodes to be the empty set
        states_present = False
        for node in tree.nodes:
            state = tree.nodedict[node].state
            if state:
                states_present = True
                tree.nodedict[node].primary_set = set([state])
                tree.nodedict[node].secondary_set = set()
        if not states_present:
            raise TreeError("There are no nodes with .state attribute in .nodedict. Can't perform parsimony analysis")

        # Tip-to-root pass: compute upper and lower sets for all unlabelled nodes
        # All leaves are assumed to be labelled (to have a non-empty .state)
        for p in tree.sorted_intnodes(deepfirst=False):
            if not tree.nodedict[p].state:
                statecount = Counter()
                for c in tree.children(p):
                    statecount.update(tree.nodedict[c].primary_set)
                statecountlist = statecount.most_common()
                topcount = statecountlist[0][1]
                maxlist = [state for state,count in statecountlist if count==topcount]
                almostmaxlist = [state for state,count in statecountlist if count==(topcount - 1)]
                tree.nodedict[p].primary_set = set(maxlist)
                tree.nodedict[p].secondary_set = set(almostmaxlist)

        return tree

    ###############################################################################################

    def _hartigan_root2tip_pass(self, tree):
        """Perform root-to-tip pass of Hartigan parsimony algorithm.
        Input: tree object with .nodedict now having .primary_set and .secondary_set
            for each node
        Output: same tree object with added attribute .optimal_set for all nodes
                The optimal_set contains all the states a node can have in a tree with MP score
                (but not all combinations of possible states will give MP score)
        """

        ndict = tree.nodedict

        # Set optimal_set of root
        ndict[tree.root].optimal_set = ndict[tree.root].primary_set

        # Root-to-tip pass: compute optimal_set for each node (labelled or not)
        for p in tree.sorted_intnodes(deepfirst=True):
            for c in tree.children(p):
                p_optimal_set = ndict[p].optimal_set
                c_primary_set = ndict[c].primary_set
                if p_optimal_set <= c_primary_set:
                    ndict[c].optimal_set = p_optimal_set
                else:
                    c_secondary_set = ndict[c].secondary_set
                    ndict[c].optimal_set = c_primary_set | (p_optimal_set & c_secondary_set)

        return tree

    ###############################################################################################

    def _set_fit_and_ambiguity_for_node(self, tree, node, parent, fitpref):
        """Set fit and ambiguity attributes for a given node."""

        ndict = tree.nodedict

        if len(ndict[node].optimal_set) == 1:
            ndict[node].wasambig = False
            ndict[node].fit = next(iter(ndict[node].optimal_set))
        else:
            ndict[node].wasambig = True
            if parent and ndict[parent].wasambig and ndict[parent].fit in ndict[node].optimal_set:
                ndict[node].fit = ndict[parent].fit
            else:
                if fitpref in ndict[node].optimal_set:
                    ndict[node].fit = fitpref
                else:
                    ndict[node].fit = random.choice(tuple(ndict[node].optimal_set))

        return tree

###################################################################################################
###################################################################################################
###################################################################################################

class TreeSet():
    """Class for storing and manipulating a number of trees, which all have the same leaves"""

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

    def nexus(self, printdist=True, printlabels=True, labelfield="label", precision=6,
                      translateblock=False, node_attributes=None, branch_attributes=None,
                      colorlist=None, colorfg="#0000FF", colorbg="#000000"):
        """Returns nexus or format tree for each tree in the set"""

        # Construct header
        stringlist = ["#NEXUS\n\n"]
        if colorlist:
            t = self[0]  # Random tree to extract leaf info
            stringlist.append(f"begin taxa\n\tdimensions ntax={len(t.leaves)};\n\ttaxlabels\n")
            for leaf in t.leaflist():
                stringlist.append(f"\t\t{leaf}")
                col = colorfg if leaf in colorlist else colorbg  # Use foreground/background colors
                stringlist.append(f"[&!color={col}]\n")
            stringlist.append(";\nend;\n\n")
        stringlist.append("begin trees;\n")

        # Add translate block if requested
        if translateblock:
            transdict = self[0].transdict()
            stringlist.append(self[0].translateblock(transdict))

        # Add newick trees with metacomments if needed
        for i, tree in enumerate(self.treelist):
            stringlist.append(f"    tree t.{i + 1} = ")
            stringlist.append(tree.newick(printdist=printdist, printlabels=printlabels, labelfield=labelfield,
                                          precision=precision, transdict=transdict if translateblock else None,
                                          node_attributes=node_attributes, branch_attributes=branch_attributes))
            stringlist.append("\n")

        # Add footer
        stringlist.append("end;\n")

        return "".join(stringlist)

    ###############################################################################################

    def newick(self, printdist=True, printlabels=True, labelfield="label", precision=6,
               transdict=None, node_attributes=None, branch_attributes=None):
        """Returns newick format tree for each tree in the set as a string"""
        stringlist = []
        for tree in self.treelist:
            stringlist.append(tree.newick(printdist=printdist, printlabels=printlabels, labelfield=labelfield,
                                          precision=precision, transdict=transdict,
                                          node_attributes=node_attributes,
                                          branch_attributes=branch_attributes))
            stringlist.append("\n")
        return "".join(stringlist)

###################################################################################################
###################################################################################################
###################################################################################################

class RootBipStruct:
    
    __slots__ = ("rootcount", "freq", "sum_frac_canon")
    
    def __init__(self, frac_to_canon=None):
        self.rootcount = 1
        self.freq = None
        self.sum_frac_canon = frac_to_canon

    def add_rootbip(self, frac_to_canon=None):
        self.rootcount += 1
        if frac_to_canon is not None:
            self.sum_frac_canon += frac_to_canon

    def merge(self, other):
        """Merges this RootBipStruct with another (for same bipartition)"""
        self.rootcount += other.rootcount
        if self.sum_frac_canon is not None:
            self.sum_frac_canon += other.sum_frac_canon

    def mean_frac_to_canon(self):
        """Return mean fraction to canonical side."""
        if self.sum_frac_canon is None:
            raise RuntimeError("Root branch-length fractions were not tracked")
        return self.sum_frac_canon / self.rootcount
                
###################################################################################################
###################################################################################################

class TreeSummary():
    """Class summarizing requested attributes (bipartitions, clades, root location, branch lengths,
       node depths, topologies) from many trees"""

    def __init__(self, trackbips=False, trackclades=False, trackroot=False, 
                       trackblen=False, trackdepth=False, trackrootblen=False,
                       tracktopo=False, track_subcladepairs=False,
                       store_trees=False):
        """TreeSummary constructor. Initializes relevant data structures"""
        self.transdict = None
        self.translateblock = None
        self.tree_count = 0
        self.tree_weight_sum = 0.0

        self.trackroot = trackroot
        self.trackbips = trackbips
        self.trackclades = trackclades
        self.trackblen = trackblen
        self.trackdepth = trackdepth
        self.trackrootblen = trackrootblen
        self.tracktopo = tracktopo
        self.track_subcladepairs = track_subcladepairs
        self.store_trees = store_trees

        self._bipartsummary = {}        # Dict: {bipartition:branchstruct with extra fields}
        self._bipartsummary_processed = False
        self._cladesummary = {}         # Dict: {clade:nodestruct with extra fields}
        self._cladesummary_processed = False
        self._rootbip_summary = {}
        self._rootbip_summary_processed = False
        self._biptoposummary = {}
        self._biptoposummary_processed = False
        self._cladetoposummary = {}
        self._cladetoposummary_processed = False

        self._sorted_biplist = None
        self._sorted_rootbips = None

    ###############################################################################################

    def __len__(self):
        return self.tree_count

    ###############################################################################################

    def online_weighted_update_mean_var(self, struct, x, w):
        """
        Weighted online update of mean + M2 (sum of weighted squared deviations).
            x: value whose mean and variance will be computed (e.g. branch length or node depth)
            w: weight

        Assumes struct has attributes: SUMW, mean, M2, n
        """

        # I am interested in being able to compute weighted mean and variance of various
        # node and branch attributes, and to do this "online" in a single pass.
        # In order to do this I follow the robust (= no underflow/overflow problems), one-pass
        # approach described in D.H.D. West, "Updating Mean and Variance Estimates: An Improved
        # Method", Communications of the ACM, 22(9), 1979.

        struct.n += 1
        delta = x - struct.mean
        struct.mean += (w / struct.SUMW) * delta
        struct.M2 += w * delta * (x - struct.mean)

    ###############################################################################################

    def finalize_online_weighted(self, struct):
        """Helper for online_weighted_update_mean_var() method:
        compute final mean, var, and sd when all values have been collected"""

        if struct.n == 0:
            msg = f"Can't compute mean and sd for item with 0 observations:\n{struct}"
            raise TreeError(msg)
        elif struct.n == 1:
            return (struct.mean, "NA", "NA")
        else:
            # Could compute variance without bias correction n/n-1 but I think mrbayes etc uses it
            var = (struct.M2 / struct.SUMW) * (struct.n / (struct.n - 1))
            sd = math.sqrt(var)
            return (struct.mean, var, sd)

    ###############################################################################################

    @property
    def bipartsummary(self):
        """Property method for lazy evaluation of freq, var, and sd for branchstructs"""
        if not self._bipartsummary_processed:
            for branchstruct in self._bipartsummary.values():
                br = branchstruct
                br.bipartition_cred = br.posterior = br.freq = br.SUMW / self.tree_weight_sum
                if self.trackblen:
                    br.length, br.length_var, br.length_sd = self.finalize_online_weighted(br)
            self._bipartsummary_processed = True
        return self._bipartsummary

    ###############################################################################################

    @property
    def cladesummary(self):
        """Property method for lazy evaluation of freq, var, and sd for nodestructs"""
        if not self._cladesummary_processed:
            for nodestruct in self._cladesummary.values():
                nd = nodestruct
                nd.clade_cred = nd.posterior = nd.freq = nd.SUMW / self.tree_weight_sum
                if self.trackdepth:
                    nd.depth, nd.depth_var, nd.depth_sd = self.finalize_online_weighted(nd)
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
                    internalbips.append((branch.posterior,bip))

            leafbips = sorted(leafbips, key=itemgetter(0))
            internalbips = sorted(internalbips, key=itemgetter(0), reverse=True)
            self._sorted_biplist = leafbips + internalbips

        return self._sorted_biplist

    ###############################################################################################

    @property
    def rootbipsummary(self):
        """Property method for lazy evaluation of freq (=rootcred) for rootbips
        Returns dict of {Bipartition: RootBipStruct}"""
        if not self._rootbip_summary_processed:
            for rootbipstruct in self._rootbip_summary.values():
                # Python note: should i divide by tree_weight_sum, not tree_count?
                rootbipstruct.freq = rootbipstruct.rootcount / self.tree_count
            self._rootbip_summary_processed = True
        return self._rootbip_summary

    ###############################################################################################

    @property
    def sorted_rootbips(self):
        """Return list of root-bipartitions (branches where root has been seen), sorted by
        occurrence (count on tree samples added to TreeSummary)"""

        if self._sorted_rootbips == None:
            self._sorted_rootbips = []
            for bip,rootbipstruct in self.rootbipsummary.items():
                self._sorted_rootbips.append((rootbipstruct._rootcount, bip, rootbipstruct))
            self._sorted_rootbips.sort(key=itemgetter(0), reverse=True)
        return self._sorted_rootbips

    ###############################################################################################

    @property
    def biptoposummary(self):
        """Property method for lazy evaluation of topostruct.posterior
        biptoposummary: dict of {biptopology: Topostruct}"""
        if not self._biptoposummary_processed:
            for topostruct in self._biptoposummary.values():
                topostruct.posterior = topostruct.weight / self.tree_weight_sum
            self._biptoposummary_processed = True

        return self._biptoposummary

    ###############################################################################################

    @property
    def cladetoposummary(self):
        """Property method for lazy evaluation of topostruct.posterior
           cladetoposummary: dict of {cladetopology: Topostruct}"""
        if not self._cladetoposummary_processed:
            for topostruct in self._cladetoposummary.values():
                topostruct.posterior = topostruct.weight / self.tree_weight_sum
            self._cladetoposummary_processed = True

        return self._cladetoposummary

    ###############################################################################################

    def add_tree(self, curtree, weight=1.0):
        """Add tree object to treesummary, update all relevant summaries"""

        # Main interface to TreeSummary.
        # Takes tree object, updates relevant measures
        # First time entered: build set of leaves for consistency checking.
        # Also compute transdict and translateblock for tree reporting (and storage?)
        if self.tree_count == 0:
            self.leaves = curtree.leaves
            self.transdict = curtree.transdict()
            self.translateblock = curtree.translateblock(self.transdict)
        elif curtree.leaves != self.leaves:
            msg = f"Leaves on tree number {self.tree_count +1} are different than previous trees"
            raise TreeError(msg)

        self.tree_count += 1
        self.tree_weight_sum += weight       # The weighted equivalent of tree_count

        if self.trackroot:
            self._add_root(curtree)

        if self.trackbips:
            self._add_bips(curtree, weight)

        if self.trackclades:
            self._add_clades(curtree, weight)

    ###############################################################################################

    def _add_root(self, curtree):
        """Helper method for add_tree: handles roots"""

        self._rootbip_summary_processed = False
        self._sorted_rootbips = None

        if self.trackrootblen:
            rootbip, frac = curtree.rootbip_with_frac()
        else:
            rootbip = curtree.rootbip()
            frac = None

        s = self._rootbip_summary.get(rootbip)
        if s is None:
            self._rootbip_summary[rootbip] = RootBipStruct(frac)
        else:
            s.add_rootbip(frac)
        
    ###############################################################################################

    def _add_clades(self, curtree, weight):
        """Helper method to add_tree: handles clades (streaming, optional topo/pairs)"""

        self._cladesummary_processed = False
        cladesummary = self._cladesummary
        trackdepth = self.trackdepth
        trackpairs = self.track_subcladepairs
        tracktopo = self.tracktopo
        online_update = self.online_weighted_update_mean_var

        # Only allocate these when needed
        node2clade = {} if trackpairs else None
        cladeset = set() if tracktopo else None

        # Stream clades: update summary + (optionally) fill node2clade + collect topology
        for node, clade, depth, nleaves in curtree.iter_cladeinfo(node2clade=node2clade):
            if cladeset is not None:
                cladeset.add(clade)

            s = cladesummary.get(clade)

            # First time clade is seen
            if s is None:
                s = Nodestruct(depth, nleaves)
                s.SUMW = weight
                if trackdepth:
                    s.n = 1
                    s.mean = depth
                    s.M2 = 0.0
                cladesummary[clade] = s

            # Clade seen before
            else:
                s.SUMW += weight
                if trackdepth:
                    online_update(s, depth, weight)

        # If requested: add subcladepairs to the GLOBAL summary
        if trackpairs:
            child_dict = curtree.child_dict
            for node in curtree.intnodes:
                try:
                    kid1, kid2 = child_dict[node].keys()
                except ValueError:
                    msg = (
                        f"Error while tracking subclade_pairs - clade has more than 2 children: "
                        f"\n\tparent-node: {node}\n\tchildren:    {curtree.children(node)}"
                        f"\nAre you sure input trees are rooted, e.g. using a clock-model?"
                    )
                    raise TreeError(msg)

                parent_clade = node2clade[node]
                parent_ns = cladesummary[parent_clade]
                if parent_ns.nleaves > 2:
                    parent_ns.add_subcladepair(node2clade[kid1],node2clade[kid2])

        # If tracking topology: update topo summary directly here (no need to return cladedict)
        if tracktopo:
            self._cladetoposummary_processed = False
            topology = frozenset(cladeset)

            ts = self._cladetoposummary.get(topology)
            if ts is None:
                ts = Topostruct()
                ts.weight = weight
                self._cladetoposummary[topology] = ts
                if self.store_trees:
                    curtree.clear_caches()
                    ts.tree = curtree
            else:
                ts.weight += weight

    ###############################################################################################

    def _add_bips(self, curtree, weight):
        """Helper method to add_tree: handles bipartitions"""

        self._bipartsummary_processed = False
        self._sorted_biplist = None
        online_update = self.online_weighted_update_mean_var
        bipartsummary = self._bipartsummary
        trackblen = self.trackblen
        tracktopo = self.tracktopo
        bipset = set() if tracktopo else None

        for bipart, branchstruct in curtree.iter_bipinfo():
            if bipset is not None:
                bipset.add(bipart)
            length = branchstruct.length
            s = bipartsummary.get(bipart)
            if s is None:
                s = branchstruct
                s.SUMW = weight
                if trackblen:
                    s.n = 1
                    s.mean = length
                    s.M2 = 0.0
                bipartsummary[bipart] = s
            else:
                s.SUMW += weight
                if trackblen:
                    online_update(s, length, weight)

        # If tracking topology, update it here
        if tracktopo:
            self._biptoposummary_processed = False
            topology = frozenset(bipset)

            ts = self._biptoposummary.get(topology)
            if ts is None:
                ts = Topostruct()
                ts.weight = weight
                self._biptoposummary[topology] = ts
                if self.store_trees:
                    curtree.clear_caches()
                    ts.tree = curtree
            else:
                ts.weight += weight

    ###############################################################################################

    def _addbiptopo(self, bipdict, curtree, weight):

        self._biptoposummary_processed = False

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

        if self.tracktopo:
            self._updatetopo(other)

    ###############################################################################################

    def _updatebip(self, other):

        # Merge "self.bipartsummary" with "other.bipartsummary"
        other_bipsum = other.bipartsummary
        self_bipsum = self.bipartsummary

        for bipart in other_bipsum:
            # If bipart already in self.bipartsummary, update fields
            if bipart in self_bipsum:
                if self.trackblen:
                    sumw1 = self_bipsum[bipart].SUMW
                    sumw2 = other_bipsum[bipart].SUMW
                    mean1 = self_bipsum[bipart].mean
                    mean2 = other_bipsum[bipart].mean
                    M21 = self_bipsum[bipart].M2
                    M22 = other_bipsum[bipart].M2
                    self_bipsum[bipart].n += other_bipsum[bipart].n
                    self_bipsum[bipart].mean = (mean1*sumw1 + mean2*sumw2)/(sumw1+sumw2)
                    self_bipsum[bipart].M2 = M21+M22+sumw1*sumw2*(mean2-mean1)*(mean2-mean1)/(sumw1+sumw2)
                self_bipsum[bipart].SUMW += other_bipsum[bipart].SUMW

            # If bipartition has never been seen before: transfer Branchstruct from other_bipsum:
            else:
                self_bipsum[bipart] = other_bipsum[bipart]

        self._bipartsummary_processed = False
        self._sorted_biplist = None

    ###############################################################################################

    def _updateclade(self, other):

        # Merge "treesummary.cladesummary" with "self.cladesummary"
        other_cladesum = other.cladesummary
        self_cladesum = self.cladesummary

        for clade in other_cladesum:
            # If bipart already in self.cladesummary, update fields
            if clade in self_cladesum:
                if self.trackdepth:
                    sumw1 = self_cladesum[clade].SUMW
                    sumw2 = other_cladesum[clade].SUMW
                    mean1 = self_cladesum[clade].mean
                    mean2 = other_cladesum[clade].mean
                    M21 = self_cladesum[clade].M2
                    M22 = other_cladesum[clade].M2
                    self_cladesum[clade].n += other_cladesum[clade].n
                    self_cladesum[clade].mean = (mean1*sumw1 + mean2*sumw2)/(sumw1+sumw2)
                    self_cladesum[clade].M2 = M21+M22+sumw1*sumw2*(mean2-mean1)*(mean2-mean1)/(sumw1+sumw2)
                self_cladesum[clade].SUMW += other_cladesum[clade].SUMW

            # If clade has never been seen before: transfer Nodestruct from other_cladesum:
            else:
                self_cladesum[clade] = other_cladesum[clade]

        self._cladesummary_processed = False
        self._sorted_biplist = None

    ###############################################################################################

    def _updateroot(self, other):

        other_rootbipsum = other._rootbip_summary
        self_rootbipsum = self._rootbip_summary
        for other_rootbip, other_s in other_rootbipsum.items():
            if other_rootbip in self_rootbipsum:
                self_rootbipsum[other_rootbip].merge(other_s)
            else:
                self_rootbipsum[other_rootbip] = other_s

    ###############################################################################################

    def _updatetopo(self, other):
        """Helper for update(). Handles biptopologies and cladetopologies"""

        self._biptoposummary_processed = False
        self._cladetoposummary_processed = False

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

    def compute_sumtree(self, treetype="con", blen="biplen", rooting=None,
                        og=None, wt_count_burnin_filename_list=None):
        """Compute and annotate summary tree: find topology, set root, set branch lengths,
           annotate branches and nodes with relevant available information (eg, sd for blen)

           Possible values (should be provided as strings):
           treetype: con (consensus tree)
                     all (consensus tree with all compatible bipartitions)
                     mcc (maximum clade credibility tree)
                     mbc (maximum bipartition credibility tree)
                     hip (HIPSTR summary tree)
                     mrhip (majority rule HIPSTR summary tree)

           blen: biplen (mean branch length for bipartitions corresponding to branches)
                 meandepth (set node depths to mean for monophyl clade, derive blens from depths)
                 cadepth (set node depths to mean common ancestor depths, derive blens from depths)
                 input (use the depths on chosen input tree, derive blens. Only valid for mcc and mbc)
                 none (all branch lengths set to 0.0)

           rooting: mid (midpoint rooting)
                    minvar (minimum variance rooting)
                    og (outgroup as string or list of strings, provide value in parameter og)
                    input (use the root on the chosen input tree - only valid for mcc and mbc trees)


           og: if rooting=og: name of outgroup taxon or list of names
        """

        # Check that all required parameters are given and consistent
        if (rooting == "og") and (not og):
            raise TreeError(f"Outgroup rooting requested, but no og parameter provided")
        if (blen == "cadepth") and (not wt_count_burnin_filename_list):
            raise TreeError("Requested cadepth but no wt_count_burnin_filename_list provided")

        # Choose summary tree topology
        if treetype == "mcc":
            sumtree = self.max_clade_cred_tree()
        elif treetype == "mbc":
            sumtree = self.max_bipart_cred_tree()
        elif treetype in ("con", "all"):
            sumtree = self.contree(allcompat=(treetype == "all"))
        elif treetype in ("hip", "mrhip"):
            sumtree = self.hipstr_tree(majrule=(treetype == "mrhip"))
        else:
            raise TreeError(f"Unknown summary tree type: {treetype}")

        # Root summary tree
        if rooting == "mid":
            sumtree.rootmid()
        elif rooting == "minvar":
            sumtree.rootminvar()
        elif rooting == "og":
            sumtree.rootout(og)
        elif rooting is None:
            pass
        else:
            raise TreeError(f"Unknown rooting method: {rooting}")

        # Set branch lengths (or depths and then branch lengths)
        if blen == "none":
            # either leave as-is or force zero
            for p in sumtree.intnodes:
                for c in sumtree.children(p):
                    sumtree.set_branch_attribute(p, c, "length", 0.0)
        elif blen == "input":
            pass  # keep whatever is already on the chosen tree
        elif blen == "meandepth":
            sumtree = self.set_mean_node_depths(sumtree)
            sumtree.set_blens_from_depths()
        elif blen == "cadepth":
            sumtree = self.set_ca_node_depths(sumtree, wt_count_burnin_filename_list)
            sumtree.set_blens_from_depths()
        elif (blen == "biplen"):
            if treetype == "mcc":
                sumtree = self.set_mean_biplen(sumtree)
        else:
            raise TreeError(f"Unknown branch-length method: {blen}")

        # Annotate tree with relevant, available attributes
        sumtree = self.annotate_sumtree(sumtree)

        # Decide what meta-comments to print
        node_attrs = set()
        branch_attrs = set()

        if self.trackclades:
            node_attrs.add("clade_cred")
        if self.trackbips:
            branch_attrs.add("bipartition_cred")

        if blen in ("meandepth", "cadepth"):
            node_attrs.update({"depth", "depth_sd", "depth_var"})

        if blen == "biplen":
            branch_attrs.add("length")
            if self.trackblen:
                branch_attrs.update({"length_sd", "length_var"})

        elif blen in ("meandepth", "cadepth"):
            branch_attrs.add("length")

        elif blen == "none":
            branch_attrs.add("length")

        elif blen == "input":
            branch_attrs.add("length")

        if self.trackroot:
            branch_attrs.add("rootcred")

        # attach for printing
        sumtree._print_node_attributes = tuple(sorted(node_attrs))
        sumtree._print_branch_attributes = tuple(sorted(branch_attrs))

        return sumtree

    ###############################################################################################

    def log_bipart_credibility(self, biptopology):
        """Compute log bipartition credibility for topology (sum of log(freq) for all branches)"""

        bipartsummary = self.bipartsummary
        logsum = 0.0
        for bipartition in biptopology:
            logsum += math.log(bipartsummary[bipartition].posterior)
        return logsum

    ###############################################################################################

    def log_clade_credibility(self, cladetopology):
        """Compute log clade credibility for topology (sum of log(freq) for all clades)"""

        cladesummary = self.cladesummary
        logsum = 0.0
        for clade in cladetopology:
            logsum += math.log(cladesummary[clade].posterior)
        return logsum

    ###############################################################################################

    def contree(self, cutoff=0.5, allcompat=False):
        """Find consensus tree built from selected bipartitions. Annotate tree with logcred"""

        if cutoff < 0.5:
            msg = "Consensus tree cutoff has to be at least 0.5"
            raise TreeError(msg)

        # Transfer biparts and branches with freq>cutoff to new bipdict, create tree
        conbipdict = {}
        i = 0
        for _, bip in self.sorted_biplist:
            i += 1
            branch = self.bipartsummary[bip]
            if branch.posterior < cutoff:
                break
            conbipdict[bip] = branch
        contree = Tree.from_biplist(conbipdict)

        # If allcompat has been requested: add remaining, compatible bipartitions to contree
        if allcompat:
            for j in range(i, len(self.sorted_biplist)):
                if contree.is_resolved():
                    break
                _,bip = self.sorted_biplist[j]
                branch = self.bipartsummary[bip]
                is_present, is_compatible, insert_tuple = contree.check_bip_compatibility(bip)
                if is_compatible and (not is_present):
                    parentnode, childnodes = insert_tuple
                    contree.insert_node(parentnode, childnodes, branch)

        logcred = self.log_bipart_credibility(contree.topology())
        contree.logcred = logcred
        contree.cred_type = "bipartition"

        return contree

    ###############################################################################################

    def max_bipart_cred_tree(self):
        """Find maximum bipartition credibility tree. Annotate tree with logcred/cred_type"""

        maxlogcred = -math.inf
        for biptopology in self.biptoposummary:
            logcred = self.log_bipart_credibility(biptopology)
            if logcred > maxlogcred:
                maxlogcred = logcred
                maxlogcredbiptopo = biptopology

        maxcredbipdict = {}
        for bipartition in maxlogcredbiptopo:
            branch = self.bipartsummary[bipartition]
            maxcredbipdict[bipartition] = branch

        # Build tree from bipartitions in new bipdict, annotate with logcred and type
        maxcredtree = Tree.from_biplist(maxcredbipdict)
        maxcredtree.logcred = maxlogcred
        maxcredtree.cred_type = "bipartition"

        return maxcredtree

    ###############################################################################################

    def max_clade_cred_tree(self):
        """Find maximum clade credibility tree. Annotate tree with logcred/cred_type"""

        maxlogcred = -math.inf
        for clade_topology in self.cladetoposummary:
            logcred = self.log_clade_credibility(clade_topology)
            if logcred > maxlogcred:
                maxlogcred = logcred
                maxlogcred_cladetopo = clade_topology

        maxcred_cladedict = {}
        for clade in maxlogcred_cladetopo:
            nodestruct = self.cladesummary[clade]
            maxcred_cladedict[clade] = nodestruct
        maxcredtree = Tree.from_cladedict(maxcred_cladedict)
        maxcredtree.logcred = maxlogcred
        maxcredtree.cred_type = "clade"

        return maxcredtree

    ###############################################################################################

    def hipstr_tree(self, majrule=False):
        """Construct HIPSTR summary tree, or mrHIPSTR tree (majrule=True)
        HIPSTR: highest independent posterior subtree reconstruction in TreeAnnotator X
        Baele et al., Bioinformatics, 2025, 41(10), https://doi.org/10.1093/bioinformatics/btaf488"""

        # For each clade, starting with the smallest:
        #     find hipstr clade_score and optimal pair of subclades
        # For clade of size 1 or 2:
        #     clade_score = log(freq)
        # For larger clades:
        #     clade_score = max[ clade_score(c1) + clade_score(c2) ]
        #                       + log(cladefreq) + majrule_reward (if freq>0.5)
        #     for all observed pairs of subclades
        cladesum = self.cladesummary
        majrule_reward = 1E10
        clades_by_size = [(nd.nleaves, clade, nd) for clade,nd in self.cladesummary.items()]
        clades_by_size = sorted(clades_by_size, key=itemgetter(0))

        for nleaves, clade, nd in clades_by_size:
            cladefreq = nd.freq
            log_cladefreq = math.log(cladefreq)

            if majrule and (cladefreq > 0.5):
                reward = majrule_reward
            else:
                reward = 0.0

            if nleaves == 1:
                nd.clade_score = log_cladefreq # Always log(1.0) = 0.0
                nd.best_pair = None
            elif nleaves == 2:
                nd.clade_score = log_cladefreq + reward
                nd.best_pair = None
            else:
                best_score = -math.inf
                best_pairs = []
                for c1,c2 in nd.subcladepairs:
                    freqsum = cladesum[c1].freq + cladesum[c2].freq # Used to break tied scores (see treeannotator)
                    clade_score = (cladesum[c1].clade_score +
                                   cladesum[c2].clade_score +
                                   log_cladefreq +
                                   reward)
                    if clade_score > best_score:
                        best_score = clade_score
                        best_pairs = [(freqsum, c1, c2)]
                    elif clade_score == best_score:
                        best_pairs.append((freqsum, c1, c2))

                if len(best_pairs) > 1:
                    best_pairs.sort(key = itemgetter(0), reverse=True)

                freqsum, best_sub1, best_sub2 = best_pairs[0]
                nd.clade_score = best_score
                nd.best_pair = (best_sub1, best_sub2)

        # Starting from root clade (all leaves): add the clades in max_subcladepair, and then their
        # max_subcladepair, etc until tree fully resolved (no deeper max_subcladepairs)
        nleaves, root_clade, root_nd = clades_by_size[-1]
        hip_clades = {root_clade:root_nd}
        stack = [root_nd]
        while stack:
            nd = stack.pop()
            if nd.best_pair:
                # Iterate over the two subclades
                for subclade in nd.best_pair:
                    nd = cladesum[subclade]
                    hip_clades[subclade] = nd
                    stack.append(nd)
        hipstr_tree = Tree.from_cladedict(hip_clades)
        hipstr_tree.clade_score = root_nd.clade_score
        hipstr_tree.logcred = self.log_clade_credibility(hipstr_tree.topology_clade)
        hipstr_tree.cred_type = "hipstr"

        return hipstr_tree

    ###############################################################################################

    def root_maxfreq(self, sumtree):
        """Uses info about root bipartitions in TreeSummary to place root on summary tree.
        Also sets tree attribute rootcred:
            probability (freq among input trees) of current location of root

        If tree has branch lengths:
        Divides length of root bipartition among two branches in accordance with average
        fraction of lengths seen for this rootbip across all trees."""

        # Starting with most frequent root location: find one that is compatible
        # Python note: should i just try number 1 on sorted list?
        if sumtree.is_bifurcation(sumtree.root):
            cur_rootbip, blenfrac_canonical_mask = sumtree.rootbip_with_frac()
        else:
            cur_rootbip = None
        for count, bip, summary_rootbipstruct in self.sorted_rootbips:
            if sumtree.bipart_is_present(bip):
                # Only reroot if tree not already rooted correctly
                if (cur_rootbip is None) or (bip != cur_rootbip):
                    parent,child = sumtree.find_bipart_nodes(bip)
                    sumtree.deroot()  # Python note: necessary?
                                      # reroot seems to assume not rooted at birfurcation
                                      # rethink reroot function and others depending on it!
                    sumtree.reroot(child, parent)

                # If branch lengths or node depths have been tracked:
                # Divide branch lengths for two rootkids according to fractions
                # seen for this rootbip across trees in ._rootbip_summary
                # if self.trackblen or self.trackdepth:
                #     kid1,kid2 = sumtree.children(sumtree.root)
                #     biplen = sumtree.nodedist(kid1, kid2)
                #     kid1_remkids = sumtree.remotechildren_dict[kid1]
                #     dist_to_kid1 = biplen * summary_rootbipstruct.avg_frac(kid1_remkids)
                #     dist_to_kid2 = biplen - dist_to_kid1
                #     sumtree.child_dict[sumtree.root][kid1].length = dist_to_kid1
                #     sumtree.child_dict[sumtree.root][kid2].length = dist_to_kid2
                # REWRITE TO USE NEW ROOTBIP!!!!!

                return sumtree

        # If we did not return by now, then bipart not in contree
        raise TreeError(f"sumtree tree not compatible with any observed root locations")

    ###############################################################################################

    def set_rootcredibility(self, sumtree, precision=6):
        """Returns sumtree with root credibilities as attributes on each branch
        rootcred = fraction of trees in input set where the root was on this branch (bipartition)
        If root was never on a branch: assign the value 0.0
        Added as attribute .rootcred to Branchstruct for branches on sumtree
        Also sets tree-attributes: 
            rootcred: the rootcred of the actually used root bipartition
            cumrootcred: sum of rootcred for all branches (bipartitions) included on tree
        """

        if not self.trackroot:
            msg = "Not possible to compute root credibilities: TreeSummary.trackroot is False"
            raise TreeError(msg)

        def set_branch_credibility(p, c):
            """Helper function to set root credibility for a branch."""
            remkids = sumtree.remotechildren_dict[c]
            bip = Bipartition.from_leafset(remkids, sumtree)
            if bip in self.rootbipsummary:
                rootcred = self.rootbipsummary[bip].freq
                sumtree.set_branch_attribute(p, c, "rootcred", rootcred)
            else:
                sumtree.set_branch_attribute(p, c, "rootcred", 0.0)

        # Find rootcred for all non-root bipartitions on sumtree
        for p in sumtree.sorted_intnodes():
            if p != sumtree.root:
                for c in sumtree.children(p):
                    set_branch_credibility(p, c)

        # Handle root bipartition separately
        p = sumtree.root
        if sumtree.is_bifurcation(p):
            rootbip = sumtree.rootbip()
            rootcred = self.rootbipsummary[rootbip].freq
            c1, c2 = sumtree.children(p)
            sumtree.set_branch_attribute(p, c1, "rootcred", rootcred)
            sumtree.set_branch_attribute(p, c2, "rootcred", rootcred)
            sumtree.rootcred = rootcred
            cumcred_correction = rootcred # Counted twice if i just sum over branches. Hackish solution...
        else:
            sumtree.rootcred = None   # sumtree is not rooted on bipartition, so no single rootcred
            for c in sumtree.children(p):
                set_branch_credibility(p, c)
            cumcred_correction = 0.0

        # Compute cumulated rootcred = sum of rootcred on all tree branches
        # Note: can be < 100% since sumtree may not contain all observed root bipartitions
        cumulated_rootcred = 0.0
        for p in sumtree.intnodes:
            for c in sumtree.children(p):
                cumulated_rootcred += sumtree.get_branch_attribute(p, c, "rootcred")
        cumulated_rootcred -= cumcred_correction  # If root at bifurcation: counted once too many...
        sumtree.cumulated_rootcred = cumulated_rootcred

        return sumtree

    ###############################################################################################

    def set_mean_node_depths(self, sumtree):
        """Set node depths on summary tree based on mean node depths for clades.
        Also set info about depth variation (sd and sem).
        All information obtained from TreeSummary.cladesummary

        NOTE 1: only meaningful if input trees are based on a clock model.
        NOTE 2: only works if all clades in sumtree have been observed at least once. The option
                will therefore not work with all rootings, and may also fail for majority rule
                consensus trees
        NOTE 3: only uses node depths from monophyletic clades (so some values may be set
                based on very few input trees)

        Node-depths of leaves will be the mean value observed across input trees.
        This only matters for leaves whose depth is being estimated (all other nodes
        will have constant depth across input trees).
        """

        try:
            for node in sumtree.nodes:
                remkids = sumtree.remotechildren_dict[node]
                clade = Clade.from_leafset(remkids, sumtree)
                nodestruct = self.cladesummary[clade]
                sumtree.set_node_attribute(node, "depth", nodestruct.depth)
                sumtree.set_node_attribute(node, "depth_sd", nodestruct.depth_sd)
        except KeyError as e:
            raise TreeError("Problem while setting mean node depths on summary tree: "
                            + "the following clade has not been observed among input trees. "
                            + "Summary tree is probably rooted differently from input trees - "
                            + "consider changing rooting option.\n"
                            + f"Unobserved clade:\n{e.args[0]}")

        return sumtree

    ###############################################################################################

    def set_ca_node_depths_orig(self, sumtree, wt_count_burnin_filename_list):
        """Set node depths on summary tree based on mean node depth of clade's MRCAs on set
        of input trees (same as "--height ca" in BEAST's treeannotator).
        This means that all input trees are used when computing mean depth for each node
        (not just the input trees where that exact monophyletic clade is present).
        Proper interpretation of node depth is then "depth of MRCA of descendant leaves".

        Node-depths of _leaves_ will be the mean value observed across input trees.
        This only matters for leaves whose depth is being estimated (all other nodes
        will have constant depth across input trees).
        """

        sumtree.clear_caches()   # Python note: ever necessary? Mostly worried about nodedepthdict

        # Find mean common ancestor depth for all internal nodes
        # Find mean node depth for leaves (for a leaf: mean depth == CA depth)
        # (most leaves will have constant depth across input trees, but if some
        # leaf dates are being estimated, then these will vary)

        # Create local bindings and precomputed list for variables used in tight loops
        online_weighted_update_mean_var = self.online_weighted_update_mean_var
        sumtree_remkids = sumtree.remotechildren_dict
        sumtree_remmask = sumtree.remotechildren_mask_dict
        sumtree_intnodes = sumtree.intnodes
        sumtree_leaves = sumtree.leaves
        leaf2mask_sum = sumtree.cached_attributes[3]

        # Make online accumulators for node stats
        acc = {}
        for node in sumtree.nodes:
            s = Nodestruct()
            s.SUMW = 0.0
            s.n = 0
            s.mean = 0.0
            s.M2 = 0.0
            acc[node] = s

        # list of (intnode, bitmask, start_leaf) queries, to iterate over for each input tree
        # start-leaf chosen randomly from remote-children of intnode
        # Python note: choice of start-leaf could be optimised to start as close to the root
        # as possible perhaps?
        intnode_queries = []
        for intnode in sumtree_intnodes:
            qmask = sumtree_remmask[intnode]
            start_leaf = next(iter(sumtree_remkids[intnode]))
            s = acc[intnode]
            intnode_queries.append((s, qmask, start_leaf))

        leaf_queries = [(acc[leaf], leaf) for leaf in sumtree_leaves]

        # Stream trees and update stats online
        for file_weight, count, burnin, filename in wt_count_burnin_filename_list:
            ntrees = count - burnin
            w_tree = file_weight / ntrees

            treefile = Treefile(filename)
            for _ in range(burnin):
                treefile.readtree(returntree=False)

            for input_tree in treefile:

                # force building remmask and nodedepthdict once per input-tree
                nodedepthdict = input_tree.nodedepthdict
                _ = input_tree.remotechildren_mask_dict      # Used by input_tree.find_mrca_mask

                # local binding for speed
                find_mrca_mask = input_tree.find_mrca_mask

                for s, query_mask, start_leaf in intnode_queries:
                    s.SUMW += w_tree
                    input_mrca = find_mrca_mask(query_mask, start_leaf)
                    depth = nodedepthdict[input_mrca]
                    online_weighted_update_mean_var(s, depth, w_tree)

                # leaves: mean depths (CA depth for singleton)
                for s, leaf in leaf_queries:
                    s.SUMW += w_tree
                    depth = nodedepthdict[leaf]
                    online_weighted_update_mean_var(s, depth, w_tree)

        # Write results onto sumtree with domain names
        for node in sumtree.nodes:
            s = acc[node]
            mean, var, sd = self.finalize_online_weighted(s)
            sumtree.set_node_attribute(node, "depth", mean)
            sumtree.set_node_attribute(node, "depth_var", var)
            sumtree.set_node_attribute(node, "depth_sd", sd)

        return sumtree

    ###############################################################################################

    def set_ca_node_depths_inline(self, sumtree, wt_count_burnin_filename_list):
        """Set node depths on summary tree based on mean node depth of clade's MRCAs on set
        of input trees (same as "--height ca" in BEAST's treeannotator).
        This means that all input trees are used when computing mean depth for each node
        (not just the input trees where that exact monophyletic clade is present).
        Proper interpretation of node depth is then "depth of MRCA of descendant leaves".

        Node-depths of _leaves_ will be the mean value observed across input trees.
        This only matters for leaves whose depth is being estimated (all other nodes
        will have constant depth across input trees).
        """

        sumtree.clear_caches()   # Python note: ever necessary? Mostly worried about nodedepthdict

        # Find mean common ancestor depth for all internal nodes
        # Find mean node depth for leaves (for a leaf: mean depth == CA depth)
        # (most leaves will have constant depth across input trees, but if some
        # leaf dates are being estimated, then these will vary)

        # Create local bindings and precomputed list for variables used in tight loops
        online_weighted_update_mean_var = self.online_weighted_update_mean_var
        sumtree_remmask = sumtree.remotechildren_mask_dict
        sumtree_intnodes = sumtree.intnodes
        sumtree_leaves = sumtree.leaves
        _, sorted_leaf_tup, _, _, _, _ = sumtree.cached_attributes

        # Helper method to pick leaf with smallest index from clade mask
        def pick_start_leaf_from_qmask(qmask):
            lsb = qmask & -qmask
            idx = lsb.bit_length() - 1
            return sorted_leaf_tup[idx]

        # Make online accumulators for node stats
        acc = {}
        for node in sumtree.nodes:
            s = Nodestruct()
            s.SUMW = 0.0
            s.n = 0
            s.mean = 0.0
            s.M2 = 0.0
            acc[node] = s

        # Group queries by start leaf
        queries_by_startleaf = defaultdict(list)
        for intnode in sumtree_intnodes:
            qmask = sumtree_remmask[intnode]
            start_leaf = pick_start_leaf_from_qmask(qmask)
            nleaves = qmask.bit_count()
            s = acc[intnode]
            queries_by_startleaf[start_leaf].append((s, qmask, nleaves))

        leaf_queries = [(acc[leaf], leaf) for leaf in sumtree_leaves]

        # Sort queries per start leaf from small to large clades
        for start_leaf, qlist in queries_by_startleaf.items():
            qlist.sort(key=lambda x: x[2])

        # Stream trees and update stats online
        # Python note: no calls to find_mrca_mask: inlined logic here instead
        # in order to enable sharing results seen while climbing ancestor path for given startleaf
        for file_weight, count, burnin, filename in wt_count_burnin_filename_list:
            ntrees = count - burnin
            w_tree = file_weight / ntrees

            treefile = Treefile(filename)
            for _ in range(burnin):
                treefile.readtree(returntree=False)

            for input_tree in treefile:
                nodedepthdict = input_tree.nodedepthdict
                remmask = input_tree.remotechildren_mask_dict
                parent_dict = input_tree.parent_dict

                # Start from parent of the leaf
                # Monotone climb: parent only moves upward
                for start_leaf, qlist in queries_by_startleaf.items():
                    parent = parent_dict[start_leaf]
                    for s, qmask, _nleaves in qlist:
                        while (remmask[parent] & qmask) != qmask:
                            parent = parent_dict[parent]
                        s.SUMW += w_tree
                        depth = nodedepthdict[parent]
                        online_weighted_update_mean_var(s, depth, w_tree)

                # leaves: mean depths (CA depth for singleton)
                for s, leaf in leaf_queries:
                    s.SUMW += w_tree
                    depth = nodedepthdict[leaf]
                    online_weighted_update_mean_var(s, depth, w_tree)

        # Write results onto sumtree with domain names
        for node in sumtree.nodes:
            s = acc[node]
            mean, var, sd = self.finalize_online_weighted(s)
            sumtree.set_node_attribute(node, "depth", mean)
            sumtree.set_node_attribute(node, "depth_var", var)
            sumtree.set_node_attribute(node, "depth_sd", sd)

        return sumtree

    ###############################################################################################

    set_ca_node_depths = set_ca_node_depths_inline

    ###############################################################################################

    def set_mean_biplen(self, sumtree):
        """Only to be used when goal is to set bipartition-based branch-length stats
        on sumtree, and those were not already set during construction.
        This happens when treetype is MCC
        
        NOTE: requires rooting and rootblen to be tracked in order to properly split branch 
              length on root bipartition"""

        if not self.trackrootblen:
            raise TreeError("Branch lengths for root-to-rootkids was not tracked - impossible to\n"
                            "use biplen on this summary tree (needed in order to properly split\n"
                            "branch length on root bipartition)")

        sumtree_root = sumtree.root
        for p in sumtree.sorted_intnodes():
            if p != sumtree_root:
                for c in sumtree.children(p):
                    remkids = sumtree.remotechildren_dict[c]
                    bip = Bipartition.from_leafset(remkids, sumtree)
                    brstruct = self.bipartsummary[bip]
                    sumtree.set_branch_attribute(p,c,"length", brstruct.length)
                    sumtree.set_branch_attribute(p,c,"length_var", brstruct.length_var)
                    sumtree.set_branch_attribute(p,c,"length_sd", brstruct.length_sd)

        # Handle root bipartition separately
        rootbip = sumtree.rootbip()
        rootbiplen = self.bipartsummary[rootbip].length
        kid1,kid2 = sumtree.child_dict[sumtree_root]
        
        summary_rootbipstruct = self._rootbip_summary[rootbip]
        frac_to_canon = summary_rootbipstruct.mean_frac_to_canon()

        kid1_mask = sumtree.remotechildren_mask_dict[kid1]
        is_kid1_canon = (kid1_mask == rootbip._mask)

        dist_to_kid1 = rootbiplen * (frac_to_canon if is_kid1_canon else (1.0 - frac_to_canon))
        dist_to_kid2 = rootbiplen - dist_to_kid1

        sumtree.child_dict[sumtree_root][kid1].length = dist_to_kid1
        sumtree.child_dict[sumtree_root][kid2].length = dist_to_kid2

        return sumtree

    ###############################################################################################

    def annotate_sumtree(self, sumtree):
        """
        Annotate sumtree nodes/branches with whatever TreeSummary tracked.
        """
        # Keep track of which attributes can be printed
        node_attrs = set()
        branch_attrs = set()

        # Node annotations from cladesummary
        if self.trackclades:
            all_leaves = sumtree.frozenset_leaves
            sorted_leafs = sumtree.sorted_leaf_tup
            leaf2index = sumtree.leaf2index

            for node in sumtree.nodes:
                clade = Clade.from_leafset(sumtree.remotechildren_dict[node], tree=sumtree)
                nd = self.cladesummary.get(clade)
                if nd is None:
                    raise TreeError("Problem while annotating summary tree:\n"
                                    + "the following clade has not been observed among input trees.\n"
                                    + "Check rooting of tree:\n"
                                    + f"{clade}")
                else:
                    # consensus clade WAS observed
                    sumtree.set_node_attribute(node, "clade_cred", getattr(nd, "clade_cred", nd.posterior))
                    node_attrs.add("clade_cred")
                    if self.trackdepth:
                        sumtree.set_node_attribute(node, "depth", nd.depth)
                        sumtree.set_node_attribute(node, "depth_sd", nd.depth_sd)
                        sumtree.set_node_attribute(node, "depth_var", nd.depth_var)
                        node_attrs.update({"depth", "depth_sd", "depth_var"})

        # Branch annotations from bipartsummary
        if self.trackbips:
            for parent in sumtree.sorted_intnodes():
                for child in sumtree.children(parent):
                    remkids = sumtree.remotechildren_dict[child]
                    bip = Bipartition.from_leafset(remkids, sumtree)
                    br = self.bipartsummary[bip]
                    sumtree.set_branch_attribute(parent, child, "bipartition_cred",
                                                 getattr(br, "bipartition_cred", br.posterior))
                    branch_attrs.add("bipartition_cred")

                    if self.trackblen:
                        sumtree.set_branch_attribute(parent, child, "length", br.length)
                        sumtree.set_branch_attribute(parent, child, "length_sd", br.length_sd)
                        sumtree.set_branch_attribute(parent, child, "length_var", br.length_var)
                        branch_attrs.update({"length", "length_sd", "length_var"})

        # Set Newick internal node labels (numbers after ')')
        # For each internal node (except root), set the label on its incoming branch.
        for parent in sumtree.sorted_intnodes():
            for child in sumtree.children(parent):
                if child in sumtree.leaves:
                    continue  # no label needed for tips

                br = sumtree.child_dict[parent][child]
                nd = sumtree.nodedict.get(child)

                # Prefer clade support if available on the node
                if nd is not None and hasattr(nd, "clade_cred"):
                    br.label = nd.clade_cred
                # Otherwise use bipartition support if available on the branch
                elif hasattr(br, "bipartition_cred"):
                    br.label = br.bipartition_cred
                else:
                    # Optional: set empty label (or leave unchanged)
                    br.label = ""

        # Root annotations
        if self.trackroot:
            self.set_rootcredibility(sumtree)
            branch_attrs.add("rootcred")

        return sumtree

    ###############################################################################################

    def set_clade_credibility(self, tree, precision=6):
        """Set clade credibility on provided target tree based on freq of clade in TreeSummary.

        NOTE: only works if all clades in tree have been observed at least once. The option
                will therefore not work with all rootings"""

        try:
            for child in (tree.intnodes - {tree.root}):
                remkids = tree.remotechildren_dict[child]
                child_clade = Clade.from_leafset(remkids, sumtree)
                clade_cred = self.cladesummary[child_clade].posterior
                parent = tree.parent(child)
                tree.setlabel(parent, child, f"{clade_cred:.{precision}g}")
        except KeyError as e:
            raise TreeError("Problem while setting clade credibililities: the following clade has not been "
                            + "observed among input trees: check rooting of tree."
                            + f"\n{e.args[0]}")

        return tree

########################################################################################
########################################################################################
###################################################################################################

class BigTreeSummary(TreeSummary):
    """DEPRECATED: thin alias left for backwards compatibility"""
    def __init__(self, store_trees=False, **kwargs):
        super().__init__(track_topology=True, store_trees=store_trees, **kwargs)


###################################################################################################
###################################################################################################
###################################################################################################

class TreefileBase():
    """Abstract base-class for representing tree file objects."""

    # Classes for specific formats inherit from this class and add extra stuff as needed.
    # NOTE: i am opening files in "read text" with encoding UTF-8. Will this work across platforms?
    # Python note: interner is used to store leaves and intnodes during reading multiple trees
    # from one treefile. If any tree is changed during reading (e.g. read and then reroot)
    # then the interned value is changed, and the next read tree will no longer match the
    # interned version (there is one more internal node for instance). This can lead to subtle
    # crashes. I think only defense is to not use interning when planning to change trees online


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
                tree = Tree.from_string(treestring)
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

        # Instantiate parser object once: this can be used repeatedly when parsing individual tree strings
        self.parser_obj = NewickStringParser(self.transdict)

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
            tree = Tree._from_string_private(self.parser_obj, treestring)
            tree.below_root = self.below_root
            return tree

###################################################################################################
###################################################################################################
###################################################################################################

class Treefile:
    """Factory for making Newick or Nexus treefile objects. Autodetects fileformat"""

    def __new__(klass, filename, interner=None):

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
            return Nexustreefile(filename, interner=interner)
        else:
            return Newicktreefile(filename, interner=interner)

###################################################################################################
###################################################################################################

class SymmetricMatrix:
    """Base class for symmetric matrices (especially distance matrices for trees or taxa).
    Only needs to set upper half, but can be accessed using any order of indices"""

    def __init__(self):
        self.dmat = None
        self.namelist = None
        self.n = None
        self.name2index = None
        self.index2name = None


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
    # pass
    # many tree samples (100,000) but small tree (41 leaves)
    tf = Nexustreefile("../sumt/tests/big_mrbayes_file.t")
    for tree in tf:
        pass
    # Fewer tree samples (4,000) but large tree (1071 leaves)
    tf = Nexustreefile("../sumt/tests/traitr_gapr_seq_realign_4.run1.t")
    for tree in tf:
        pass

###################################################################################################

if __name__ == "__main__":
    main()
    # import cProfile
    # cProfile.run('main()', 'tmp/profile.pstats')

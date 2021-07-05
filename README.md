# phylotreelib: python module for manipulating and analyzing phylogenetic trees

Using classes and methods in phylotreelib.py it is possible to read treefiles in either
NEXUS or Newick format, and to analyze and manipulate the trees in various ways.

![](https://github.com/agormp/phylotreelib/raw/master/phylogenetictree.png?raw=true)

## Availability

The phylotreelib.py module is available on GitHub: https://github.com/agormp/phylotreelib and on PyPI: https://pypi.org/project/phylotreelib/

## Installation

```
pip3 install phylotreelib
```

## Quick start example

Here is a script that will import phylotreelib, read a NEXUS file, perform minimum variance rooting, find the node ID for the new rootnode,
and finally print out the root-to-tip distance (measured along the branches) for all tips in the tree:

```python
import phylotreelib as pt
treefile = pt.Nexustreefile("mytreefile.nexus")
mytree = next(treefile)
mytree.rootminvar()
rootnode = mytree.root
for tip in mytree.leaves:
	dist = mytree.nodedist(rootnode, tip)
	print(dist)
```

## Using phylotreelib

The main idea is to construct a treefile object from text (e.g., obtained by
reading a NEXUS or Newick treefile). Treefile objects contain tree objects, and
can be iterated over. Tree objects have a number of methods that can be used to
analyze or alter the tree in question.

### Constructing treefile objects

Typically, treefile objects are constructed by reading a textfile containing a
NEXUS or Newick format description of a tree. This is done as follows:

Newick files:
```
treefile = phylotreelib.Newicktreefile(filename)
```

NEXUS files:
```
treefile = phylotreelib.Nexustreefile(filename)
```

### Constructing tree objects

Treefile objects contain tree objects. Tree objects are retrieved from treefile
objects by iteration:

Doing something to all trees in a treefile object:
```
for tree in treefile:
	<do something with tree>
```

Getting a single tree from a treefile object:
```
tree = next(treefile)
```

A tree object can also be constructed directly from a string (where the string is a Newick formatted tree):
```
tree = phylotreelib.Tree.from_string(mystring)
```

Tree objects consist of external nodes (leafs), which are identified by strings
(e.g. "Chimpanzee"), and internal nodes, which are identified by integers
(e.g., 5). Branches between nodes may have a label (string) and/or a branch
length (float) associated with them.

A textual representation of a tree object can be obtained using: "print(mytree)" (where mytree is the name of the tree object). The resulting
output is a child-list representation of the tree followed by an alphabetical list of leafs, along these lines:

```
>>> print(tree)
|------------------------------------------------------------------|
|  Node  |              Child               |  Distance  |  Label  |
|------------------------------------------------------------------|
|     0  |  tetM_AF333235_Clostridium_diff  |  0.045404  |         |
|     0  |                         9830414  |  0.003285  |         |
|     0  |                         9830359  |  0.001099  |         |
|     0  |                               1  |  0.018480  |   1.00  |
|     1  |                              24  |  0.007180  |   1.00  |
|     1  |                               2  |  0.026523  |   1.00  |
|    24  |                         9830470  |  0.001096  |         |
|    24  |                         9830409  |  0.000563  |         |
|     2  |                               3  |  0.016300  |   0.77  |
|     2  |  tetO_AY190525_Campylobacter_je  |  0.565353  |         |
|     3  |                              17  |  0.007065  |   0.82  |
|     3  |                               4  |  0.006912  |   1.00  |
|    17  |                              18  |  0.004701  |   1.00  |
|    17  |                              19  |  0.001937  |   0.87  |
|     4  |                               5  |  0.003834  |   1.00  |
|     4  |                               6  |  0.007932  |   1.00  |
|    18  |  tetM_X90939_Streptococcus_pneu  |  0.000804  |         |
|    18  |   tetM_AF376746rc_Streptococcus  |  0.001385  |         |
|    19  |                              20  |  0.001780  |   0.98  |
|    19  |                              21  |  0.007687  |   0.97  |
|     5  |                           20028  |  0.001084  |         |
|     5  |                           18854  |  0.000562  |         |
|     6  |  tetM_U58986_Gardnerella_vagina  |  0.007218  |         |
|     6  |                               7  |  0.010529  |   1.00  |
|    20  |  tetM_AJ580978_Streptococcus_mi  |  0.004311  |         |
|    20  |  tetM_AJ580977_Streptococcus_mi  |  0.002311  |         |
|    21  |  tetM_L12242_Neisseria_gonorrho  |  0.001760  |         |
|    21  |                              22  |  0.011284  |   1.00  |
|     7  |  tetM_L12241_Neisseria_gonorrho  |  0.001992  |         |
|     7  |                               8  |  0.008244  |   1.00  |
|     8  |                               9  |  0.006822  |   1.00  |
|     8  |  tetM_U58985_Gardnerella_vagina  |  0.000684  |         |
|    22  |  tetM_X04388_Enterococcus_faeca  |  0.000610  |         |
|    22  |                              23  |  0.002168  |   1.00  |
|     9  |                              10  |  0.002313  |   1.00  |
|     9  |                              12  |  0.001110  |   0.83  |
|    23  |                           20074  |  0.000577  |         |
|    23  |  tetM_AB054984_Clostridium_sept  |  0.002784  |         |
|    10  |                           20032  |  0.000537  |         |
|    10  |                              11  |  0.002783  |   1.00  |
|    10  |  tetM_X92947_Enterococcus_faeca  |  0.001092  |         |
|    10  |                         9830090  |  0.000543  |         |
|    10  |                         9830498  |  0.000536  |         |
|    10  |  tetM_M85225_Enterococcus_faeca  |  0.000536  |         |
|    10  |  tetM_U09422_Enterococcus_faeca  |  0.000541  |         |
|    10  |                         9830479  |  0.000547  |         |
|    10  |                           18836  |  0.000542  |         |
|    10  |                           15109  |  0.000540  |         |
|    10  |                         9830457  |  0.000536  |         |
|    10  |                         9830491  |  0.000537  |         |
|    12  |                              16  |  0.001123  |   0.99  |
|    12  |                              13  |  0.005133  |   1.00  |
|    12  |                              14  |  0.002217  |   1.00  |
|    16  |  tetM_ABO39845_Erysipelothrix_r  |  0.000567  |         |
|    16  |  tetM_AF329848_Clostridium_perf  |  0.000547  |         |
|    11  |  tetM_X56353_Enterococcus_faeca  |  0.000559  |         |
|    11  |  tetM_AJ580976_Streptococcus_or  |  0.004525  |         |
|    13  |  tetM_U08812_Ureaplasma_urealyt  |  0.004481  |         |
|    13  |  tetM_X75073_Neisseria_meningit  |  0.001664  |         |
|    13  |  tetM_AOE14233_Streptococcus_ag  |  0.000558  |         |
|    13  |  tetM_AF440277_Lactobacillus_pl  |  0.001113  |         |
|    14  |   tetM_AF491293_Bacillus_cereus  |  0.000589  |         |
|    14  |                              15  |  0.022849  |   1.00  |
|    15  |  tetM_M21136_Staphylococcus_aur  |  0.001704  |         |
|    15  |  tetM_APO03359rc_Staphylococcus  |  0.000573  |         |
|------------------------------------------------------------------|

41 Leafs:
------------------------------
15109
18836
18854
20028
20032
20074
9830090
9830359
9830409
9830414
9830457
9830470
9830479
9830491
9830498
tetM_AB054984_Clostridium_sept
tetM_ABO39845_Erysipelothrix_r
tetM_AF329848_Clostridium_perf
tetM_AF333235_Clostridium_diff
tetM_AF376746rc_Streptococcus
tetM_AF440277_Lactobacillus_pl
tetM_AF491293_Bacillus_cereus
tetM_AJ580976_Streptococcus_or
tetM_AJ580977_Streptococcus_mi
tetM_AJ580978_Streptococcus_mi
tetM_AOE14233_Streptococcus_ag
tetM_APO03359rc_Staphylococcus
tetM_L12241_Neisseria_gonorrho
tetM_L12242_Neisseria_gonorrho
tetM_M21136_Staphylococcus_aur
tetM_M85225_Enterococcus_faeca
tetM_U08812_Ureaplasma_urealyt
tetM_U09422_Enterococcus_faeca
tetM_U58985_Gardnerella_vagina
tetM_U58986_Gardnerella_vagina
tetM_X04388_Enterococcus_faeca
tetM_X56353_Enterococcus_faeca
tetM_X75073_Neisseria_meningit
tetM_X90939_Streptococcus_pneu
tetM_X92947_Enterococcus_faeca
tetO_AY190525_Campylobacter_je
```

### Attributes on Tree objects

Tree objects have a number of attributes that contain information regarding the tree.

One example of using a tree attribute is (assuming we have a Tree object named "tree"):

```
taxa = tree.leaves
```

List of  useful Tree object attributes:

*	`.leaves` 	: Set of leaf names
*	`.intnodes` 	: Set of internal node IDs
*	`.nodes`		: Set of all nodes (= leafs + internal nodes)
*	`.root`		: ID for root node (usually 0, but may change if re-rooted)


### Methods for analyzing and altering tree objects.

Tree objects also have a number of methods that can be used to analyze and alter them. In
addition to what's mentioned below, phylotreelib also contains additional methods, but
these are mostly for internal use in phylotreelib (not for application programming),
and are subject to change if I come up with better implementations.

One example of using a tree object method is:

```
childnodes = tree.children(7)
```

which returns a set containing the node-IDs for the immediate descendants of node 7 (the nodes directly connected to node 7).

A full list of classes and methods in phylotreelib is at the end of this README

### Exceptions.

phylotreelib has its own error class ("TreeError"), which may be handy for catching
tree-related errors in your own program and dealing with it intelligently:

Example usage:

```
try:
    tree.rootout(outgroup)
except phylotreelib.TreeError as err:
	print("This error occurred: {}".format(err) )
```

Example usage 2:

```
if nodename not in tree.nodes:
	raise TreeError("Tree contains no leafs named {}".format(nodename))
```


### List of  methods in the Tree class


```
class Tree(builtins.object)
 |  Class representing basic phylogenetic tree object.
 |
 |  Methods defined here:
 |
 |  __eq__(self, other)
 |      Implements equality testing for Tree objects
 |
 |  __hash__(self)
 |      Implements hashing for Tree objects, so they can be used as keys in dicts
 |
 |  __init__(self)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |
 |  __iter__(self)
 |      Returns iterator object for Tree object. Yields subtrees with extra .basalbranch attribute (= basal branch struct)
 |
 |  __str__(self)
 |      Prints a table of parent-child relationships in the tree including branch lengths and labels
 |
 |  add_branch(self, bipart, blen=0.0, label='')
 |      Adds branch represented by bipartition to unresolved tree.
 |
 |  average_ancdist(self, leaflist, return_median=False)
 |      Return average distance from leaves to their MRCA, measured along tree. Median can be requested
 |
 |  average_pairdist(self, leaflist, return_median=False)
 |      Return average pairwise distance between leaves in leaflist, measured along tree. Median can be requested
 |
 |  bipdict(self)
 |      Returns tree in the form of a "bipartition dictionary"
 |
 |  build_dist_dict(self)
 |      Construct dictionary keeping track of all pairwise distances between nodes
 |
 |  build_parent_dict(self)
 |      Forces construction of parent_dict enabling faster lookups (avoid function, use dict directly)
 |
 |  build_path_dict(self)
 |      Construct dictionary keeping track of all pairwise paths between nodes
 |
 |  children(self, parent)
 |      Returns set containing parent's immediate descendants
 |
 |  cladegrep(self, pattern, minsize=2)
 |      Finds clades (monophyletic groups) where all leaves contain specified pattern
 |
 |  cluster_cut(self, cutoff)
 |      Divides tree into clusters by conceptually cutting across tree at "cutoff" distance from root.
 |      Returns list containing sets of leafnames
 |
 |  cluster_n(self, nclust, return_as_basenodes=False)
 |      Divides tree into 'nclust' clusters based on distance from root.
 |      Returns list containing sets of leafnames unless requested to return basenodes
 |
 |  deroot(self)
 |      If root is at bifurcation: remove root node, connect adjacent nodes
 |
 |  diameter(self, return_leaves=False)
 |      Return diameter: longest leaf-leaf distance along tree. If return_leaves is True: Return tuple with (maxdist, Leaf1, Leaf2)
 |
 |  figtree(self, printdist=True, printlabels=True, precision=6, colorlist=None, color='0000FF')
 |      Returns figtree format tree as a string. Rudimentary - mostly for coloring leaves. Default color=blue
 |
 |  findMRCA(self, leafset)
 |      Finds Most Recent Common Ancestor for the provided set of leaves
 |
 |  find_central_leaf(self, leaflist)
 |      Finds central leaf for the provided list of leaves.
 |      Defined as having approximately equal distance to the two farthest leaves in leaflist
 |
 |  find_common_leaf(self, leaflist)
 |      Finds common leaf for the provided list of leaves.
 |      Defined as having the smallest average distance to remaining leaves (= many close neighbors).
 |
 |  find_most_distant(self, node1, nodeset)
 |      Finds node in nodeset that is most distant from node1
 |
 |  findbasenode(self, leafset)
 |      Finds node that is at the base of all leaves in leafset.
 |
 |  getlabel(self, node1, node2)
 |      Gets label on branch connecting node1 and node2
 |
 |  graft(self, other, node1, node2=None, blen1=0, blen2=0, graftlabel=None)
 |      Graft other tree to self. Tree2 (other) intnodes renamed if names clash with those in tree1.
 |      node1: node in tree1 (self) below which tree2 (other) will be grafted. Must be specified, cannot be root1
 |      node2: node in tree2 (other) below which tree2 will be attached (defaul is root of tree2)
 |      blen1: length of branch added to tree1 below graftpoint (lower of two newly created branches)
 |      blen2: length of branch above graft point and below tree2 (upper of two newly created branches)
 |      graftlabel: prepend value of "label" to leaf names on t2 (e.g: "graft_s1")
 |
 |  insert_node(self, parent, childnodes, branchlength=0, lab='')
 |      Inserts an extra node between parent and children listed in childnodes list.
 |
 |      Length of inserted branch is 'branchlength' and defaults to zero. The node number
 |      of the new node is returned
 |
 |  is_compatible_with(self, bipart)
 |      Checks whether a given bipartition is compatible with the tree
 |
 |  is_resolved(self)
 |      Checks whether tree is fully resolved (no polytomies)
 |
 |  leaflist(self)
 |      Returns list of leaf names
 |
 |  length(self)
 |      Returns tree length (sum of all branch lengths)
 |
 |  nameprune(self, sep='_', keep_pattern=None)
 |      Prune leaves based on name redundancy: Find subtrees where all leaves have same start of name (up to first "_")
 |
 |  newick(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6)
 |      Returns Newick format tree string representation of tree object
 |
 |  nexus(self, printdist=True, printlabels=True, precision=6)
 |      Returns nexus format tree as a string
 |
 |  nodedepth(self, node)
 |      Returns depth of node: distance from rightmost leaf-level to node (i.e., depth of rightmost leaf = 0.0)
 |
 |  nodedist(self, node1, node2=None)
 |      Returns distance between node1 and node2 along tree (patristic distance)
 |
 |  nodepath(self, node1, node2)
 |      Returns path between node1 and node2 along tree.
 |
 |  nodepath_fromdict(self, node1, node2)
 |      Returns path between node1 and node2 along tree, from preconstructed path_dict
 |
 |  numberprune(self, nkeep, keeplist=None, keep_common_leaves=False, keep_most_distant=False, return_leaves=False, enforceN=False)
 |      Prune tree so 'nkeep' leaves remain. Leaves are chosen to be approximately evenly spaced over tree.
 |      "keeplist" can be used to specify leaves that _must_ be retained.
 |      'keep_common_leaves' requests preferential retainment of leaves with many neighbors
 |      (default is to keep leaves that are as equally spaced as possible)
 |      'keep_most_distant' requests that the two most distant leaves in tree (which spread out the diameter) should be kept
 |      'return_leaves': return selected leaves, but do not actually prune tree
 |      'enforceN' enforce exactly N leaves in pruned tree (normally leaves in includelist and most distant are additional to N)
 |
 |  parent(self, node)
 |      Returns parent of node
 |
 |  prune_maxlen(self, nkeep, return_leaves=False)
 |      Prune tree so the remaining nkeep leaves spread out maximal percentage of branch length (max phylogenetic diversity)
 |
 |  remote_children(self, parent)
 |      Returns set containing all leaves that are descendants of parent
 |
 |  remove_branch(self, node1, node2)
 |      Removes branch connecting node1 and node2 (thereby creating polytomy)
 |
 |  remove_leaf(self, leaf)
 |      Removes named leaf from tree, cleans up so remaining tree structure is sane
 |
 |  remove_leaves(self, leaflist)
 |
 |  rename_intnode(self, oldnum, newnum)
 |      Changes number of one internal node
 |
 |  rename_leaf(self, oldname, newname, fixdups=False)
 |      Changes name of one leaf. Automatically fixes duplicates if requested
 |
 |  reroot(self, node1, node2, polytomy=False, node1dist=0.0)
 |      Places new root on branch between node1 and node2, node1dist from node1
 |
 |  resolve(self)
 |      Randomly resolves multifurcating tree by by adding zero-length internal branches.
 |
 |  rootmid(self)
 |      Performs midpoint rooting of tree
 |
 |  rootout(self, outgroup, polytomy=False)
 |      Roots tree on outgroup
 |
 |  setlabel(self, node1, node2, label)
 |      Sets label on branch connecting node1 and node2
 |
 |  setlength(self, node1, node2, length)
 |      Sets length of branch connecting node1 and node2
 |
 |  shuffle_leaf_names(self)
 |      Shuffles the names of all leaves
 |
 |  sorted_intnodes(self, deepfirst=True)
 |      Returns sorted intnode list for breadth-first traversal of tree. Default is to place deep nodes first
 |
 |  spr(self, subtree_node, regraft_node)
 |      Subtree Pruning and Regrafting.
 |      subtree_node: basenode of subtree that will be pruned.
 |      regraft_node: node in tree below which subtree will be grafted. Must be specified, cannot be root
 |
 |  subtree(self, basenode, return_basalbranch=False)
 |      Returns subtree rooted at basenode as Tree object. Note: rooting matters! Note 2: basenode may be leaf!
 |
 |  topology(self)
 |      Returns set of sets of sets representation of topology ("naked bipartitiondictionary")
 |
 |  transname(self, namefile)
 |      Translate all leaf names using oldname/newname pairs in namefile
 |
 |  treedist(self, other, normalise=True, verbose=False)
 |      Compute symmetric tree distance (Robinson Foulds) between self and other tree. Normalised measure returned by default
 |
 |  treesim(self, other, verbose=False)
 |      Compute normalised symmetric similarity between self and other tree
 |
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |
 |  from_string(orig_treestring, transdict=None) from builtins.type
 |      Constructor 1: Converts tree (string) in Newick format to internal representation
 |
 |  from_biplist(biplist) from builtins.type
 |      Constructor 2: constructs Tree object from bipartition list
 |
 |  from_leaves(leaflist) from builtins.type
 |      Constructor 3: construct star-tree object from list of leaves
 |
 |  randtree(leaflist=None, ntips=None, randomlen=False, name_prefix='s') from builtins.type
 |      Constructor 4: Construct tree with random topology. Either list of leaf names OR number of tips must be specified
 |
 |  ----------------------------------------------------------------------

```


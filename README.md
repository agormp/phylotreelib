# phylotreelib: python library for analyzing and manipulating phylogenetic trees

[![PyPI downloads](https://static.pepy.tech/personalized-badge/phylotreelib?period=total&units=none&left_color=black&right_color=blue&left_text=downloads&service=github)](https://pepy.tech/project/phylotreelib)
![](https://img.shields.io/badge/version-1.6.0-blue)


Using classes and methods in phylotreelib.py it is possible to read and write treefiles and to analyze and manipulate the trees in various ways.

![](https://github.com/agormp/phylotreelib/raw/master/phylogenetictree.png?raw=true)

## Availability

The phylotreelib.py module is available on GitHub: https://github.com/agormp/phylotreelib and on PyPI: https://pypi.org/project/phylotreelib/

## Installation

```
python3 -m pip install phylotreelib
```

## Highlights

* Methods for reading and writing Nexus and Newick format tree files
* Trees can be interrogated to find e.g. the parents or children of specific nodes, to find the most recent common ancestor (MRCA) of a set of leaves, to find the N leaves that are closest to a specified leaf, to find the longest path on the tree, etc.
* Possible to iterate over all subtrees in Tree object
* Trees can be rooted using either midpoint, outgroup, or minimum-variance rooting
* Trees can be pruned to a smaller, specified number of representative leaves
* Nodes can be inserted or deleted while keeping the rest of the tree structure sane, leaves can be renamed, branch lengths can be set
* Trees can be grafted to each other, or a subtree can be moved using subtree pruning and regrafting (SPR)
* Method for computing the distance between two trees (Robinson-Foulds, and derived normalised measures; returns information about number of unique and shared bipartitions)
* The Distmatrix class contains methods for building trees from distance matrices
* Methods for computing consensus trees from sets of input trees, which also yield information on the frequencies of clades and topologies, and on the distribution of branch lengths for bipartitions
* Library has been optimized for high speed and low memory consumption

## Quick start usage examples

### Read tree from file, perform minimum-variance rooting, find root-to-tip distances
The code below will import phylotreelib, open a NEXUS file, read one Tree object from the file, perform minimum-variance rooting, find the node ID for the new rootnode,
and finally print out the name and root-to-tip distance (measured along the branches) for all leaves in the tree:

```python
import phylotreelib as pt
treefile = pt.Nexustreefile("mytreefile.nexus")
mytree = treefile.readtree()
mytree.rootminvar()
rootnode = mytree.root
for tip in mytree.leaves:
	dist = mytree.nodedist(rootnode, tip)
	print("{:<10s}\t{:.2f}".format(tip,dist))
```

Output:

```
nitrificans	1879.84
Is79A3    	1878.95
GWW4      	1877.47
.
.
.
A2        	1879.84
communis  	1878.95
```

-----

### Construct tree from Newick string, print tabular representation
The code below constructs a Tree object from a Newick formatted string and then prints the string representation of the tree (using the Tree object's \_\_str\_\_() method).

```python
import phylotreelib as pt
mytree = pt.Tree.from_string("(Gorilla:3, (Human:2, (Chimpanzee:1, Bonobo:1):1):1);")
print(mytree)
```

Output:

```
|----------------------------------------------|
|  Node  |    Child     |  Distance  |  Label  |
|----------------------------------------------|
|     0  |           1  |         1  |         |
|     0  |     Gorilla  |         3  |         |
|     1  |           2  |         1  |         |
|     1  |       Human  |         2  |         |
|     2  |      Bonobo  |         1  |         |
|     2  |  Chimpanzee  |         1  |         |
|----------------------------------------------|

4 Leafs:
----------
Bonobo
Chimpanzee
Gorilla
Human
```

-----

### Read multiple trees from Nexus file, construct consensus tree, write in Newick format
The code below opens a Nexus-formatted file with multiple trees, constructs a Treesummary object, and then extracts all Tree objects from the file by iterating over the file while adding the trees to the Treesummary object. Then a majority rule consensus tree is computed from the Treesummary object, the tree is midpoint rooted, and the resulting tree is finally written in Newick format to the output file "contree.newick"

```python
import phylotreelib as pt
beastfile = pt.Nexustreefile("BEAST_samples.trees")
treesummary = pt.TreeSummary()
for tree in beastfile:
    treesummary.add_tree(tree)
consensus_tree = treesummary.contree()
consensus_tree.rootmid()
with open("contree.newick", "w") as outfile:
    outfile.write(consensus_tree.newick())
```

-----

### Read tree from file, find close neighbours of specified leaf
The code below opens a Newick file, reads one Tree object from the file, and then finds the 5 leaves that are closest (measured along the branches) to the leaf labeled "nitrificans".

```python
import phylotreelib as pt
treefile = pt.Newicktreefile("Comammox.newick")
tree = treefile.readtree()
print(tree.nearest_n_leaves("nitrificans", 5))
```

Output:

```python
{'A2', 'AAUMBR1', 'CG24B', 'inopinata', 'nitrosa'}
```

-----

### Read DNA alignment, construct Neighbor Joining tree
The code below opens a fasta file containing a set of aligned DNA sequences and reads the aligned sequences (using classes and methods from the [sequencelib](https://github.com/agormp/sequencelib) library), constructs a nested dictionary containing all pairwise sequence distances, constructs a Distmatrix object from this dictionary, and computes a neighbor joining tree from the distance matrix.

```python
import phylotreelib as pt
import sequencelib as seqlib
seqfile = seqlib.Seqfile("myalignment.fasta")
seqs = seqfile.read_alignment()
distdict = seqs.distdict()
dmat = pt.Distmatrix.from_distdict(distdict)
mytree = dmat.nj()
```

## Using phylotreelib

Typically, phylotreelib will be used for analysing (or manipulating) one or more trees that have been read from a textfile in Newick or Nexus format. Reading a tree from file will return a Tree object, which has methods for interrogating or altering itself (e.g. `mytree.rootmid()` will midpoint root the Tree object `mytree`).

### Opening a treefile

To open a Newick format file:

```
treefile = phylotreelib.Newicktreefile(filename)
```

To open a Nexus format file:

```
treefile = phylotreelib.Nexustreefile(filename)
```

These commands will return a file object with methods for reading the contents of the file.

### Reading one or more trees from a treefile

To read one tree from an opened treefile:

```
tree = treefile.readtree()
```

This returns a Tree object, which has methods for analysing and manipulating the tree (itself). By calling readtree() repeatedly, you can read additional trees from the file. The `readtree()` method returns `None` when all trees have been read.

The `readtrees()` method returns all the trees in the file in the form of a Treeset object. Treeset objects contains a list of Tree objects and has methods for rooting and outputting all trees in the collection. Iterating over a Treeset object returns Tree objects.

```
treeset = treefile.readtrees()
```

It is also possible to retrieve all the trees from an open treefile one at a time by iterating directly over the file (useful for minimizing memory consumptiom when handling files with many trees):

```
for tree in treefile:
	<do something with tree>
```

### Constructing Tree objects directly

Instead of reading a tree from a file, you can also construct Tree objects using one of the 3 alternative constructors in the Tree class.

#### Constructing Tree object from string

Tree objects can be constructed directly from a string (where the string is a Newick formatted tree):

```
tree = phylotreelib.Tree.from_string(mystring)
```

#### Constructing Tree object with star-tree topology from list of leaf names

Tree objects (with a star topology) can be constructed from a list of leaf names:

```
tree = phylotreelib.Tree.from_leaves(leaflist)
```

#### Constructing Tree object with random topology and branch lengths

It is possible to construct Tree objects with random tree topology using the randtree constructor:

```
tree = phylotreelib.Tree.randtree(leaflist=None, ntips=None, randomlen=False, name_prefix="s"):
```

Either a list of names (`leaflist`) or the number of tips (`ntips`) can be specified as a way of setting the size of the tree. If the function argument `randomlen` is True then branches will get random lengths drawn from a lognormal distribution.

### Structure of Tree objects

Tree objects consist of external nodes (leaves), which are identified by strings
(e.g. "Chimpanzee"), and internal nodes, which are identified by integers
(e.g., 5). Branches between nodes may have a label (string) and/or a branch
length (float) associated with them.

A textual representation of a Tree object can be obtained using: "print(mytree)" (where mytree is the name of the Tree object). The resulting output is a child-list representation of the tree followed by an alphabetical list of leafs, along these lines:

```
>>> print(mytree)
|-------------------------------------------------|
|  Node  |     Child     |   Distance   |  Label  |
|-------------------------------------------------|
|     0  |      lmo0024  |   0.0246752  |         |
|     0  |            1  |    0.972917  |   1.00  |
|     0  |      lin0023  |   0.0627209  |         |
|     1  |     SMU_1957  |     1.02145  |         |
|     1  |            2  |    0.219277  |   0.83  |
|     2  |           17  |    0.234726  |   0.98  |
|     2  |            3  |    0.089033  |   0.81  |
|    17  |      CPE0323  |    0.595102  |         |
|    17  |      CPE2630  |    0.889163  |         |
|     3  |            4  |    0.246187  |   0.86  |
|     3  |           14  |    0.284347  |   1.00  |
|     4  |            5  |    0.124892  |   0.99  |
|     4  |            7  |    0.285426  |   0.85  |
|    14  |           16  |    0.862141  |   1.00  |
|    14  |           15  |    0.458832  |   1.00  |
|     5  |  Bsu_COG3716  |    0.574712  |         |
|     5  |            6  |    0.276761  |   1.00  |
|     7  |            8  |   0.0462592  |   0.55  |
|     7  |     SMU_1879  |     0.22478  |         |
|     7  |         Lme1  |     0.29347  |         |
|    16  |          Sag  |    0.300177  |         |
|    16  |      CPE1463  |    0.178835  |         |
|    15  |      lmo2000  |     0.10323  |         |
|    15  |      lin2108  |   0.0582161  |         |
|     6  |      lmo0781  |   0.0465803  |         |
|     6  |      lin0774  |   0.0163573  |         |
|     8  |            9  |   0.0570041  |   0.82  |
|     8  |         Lls_  |    0.318751  |         |
|     9  |           10  |   0.0290423  |   0.70  |
|     9  |       LJ0505  |    0.256245  |         |
|    10  |         Lme2  |     1.73392  |         |
|    10  |           11  |   0.0534659  |   0.54  |
|    10  |          Lpl  |    0.183936  |         |
|    10  |           13  |    0.107465  |   0.99  |
|    11  |       EF0022  |   0.0951663  |         |
|    11  |           12  |    0.122973  |   0.92  |
|    13  |      CPE0823  |    0.176274  |         |
|    13  |          CAC  |    0.113429  |         |
|    12  |      lin0145  |  0.00949095  |         |
|    12  |      lmo0098  |   0.0119641  |         |
|-------------------------------------------------|

23 Leafs:
-----------
Bsu_COG3716
CAC
CPE0323
CPE0823
CPE1463
CPE2630
EF0022
LJ0505
Lls_
Lme1
Lme2
Lpl
SMU_1879
SMU_1957
Sag
lin0023
lin0145
lin0774
lin2108
lmo0024
lmo0098
lmo0781
lmo2000
```

### Attributes on Tree objects

Tree objects have a number of attributes that contain information regarding the tree.

One example of using a tree attribute is (assuming we have a Tree object named "tree"):

```
taxa = tree.leaves
```

List of useful Tree object attributes:

*	`.leaves` 	: Set of leaf names
*	`.intnodes` : Set of internal node IDs
*	`.nodes`	: Set of all nodes (= leafs + internal nodes)
*	`.root`		: ID for root node (usually 0, but may change if re-rooted)


### Methods for analyzing and altering Tree objects.

Tree objects also have a number of methods that can be used to analyze and alter them.

One example of using a Tree object method is:

```
childnodes = tree.children(7)
```

which returns a set containing the node-IDs for the immediate descendants of node 7 (the nodes directly connected to node 7).

A full list of classes and methods in phylotreelib is at the end of this README

### Exceptions.

phylotreelib has its own error class ("TreeError"), which can be used for raising and catching
tree-related errors in your own program and dealing with errors intelligently:

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

### Most important classes in phylotreelib

* **Tree**: Class representing basic phylogenetic Tree object. Methods for analyzing and manipulating tree (rooting, finding MRCA, finding children of nodes, ...)
* **TreeSet**: Class for storing and manipulating a number of trees.
* **TreeSummary**: Class summarizing bipartitions and branch lengths (but not topologies) from many trees. Method to compute consensus tree. Can return information about clade frequencies and average branch lengths (where a branch is defined as being between the two sets of leaves in a bipartition, regardless of the topologies of the subtrees).
* **BigTreeSummary**: Class summarizing bipartitions, branch lengths, and topologies from many trees (where topology is defined as the set of bipartitions corresponding to the unrooted tree). Does everything TreeSummary does and also keeps track of topologies.
* **Nexustreefile**: Class representing Nexus tree file. Iteration returns Tree objects.
* **Newicktreefile**: Class representing Newick tree file. Iteration returns Tree objects.
* **Distmatrix**: Class representing distance matrix for set of taxa. Can be constructed from distance data or directly from alignment of sequences (requires seqlib.py). Methods for computing neighbor-joining and upgma trees."""


### Help for all classes and methods in phylotreelib

(Output from pydoc on module, classes listed in alphabetical order)

```
Help on module phylotreelib:

NAME
    phylotreelib - Classes and methods for analyzing, manipulating, and building phylogenetic trees

CLASSES
    builtins.Exception(builtins.BaseException)
        TreeError
    builtins.object
        Branchstruct
        Distmatrix
        Globals
        Topostruct
        Tree
        TreeSet
        TreeSummary
            BigTreeSummary
        Treefile
            Newicktreefile
            Nexustreefile

    class BigTreeSummary(TreeSummary)
     |  BigTreeSummary(include_zeroterms=False, outgroup=None, rootmid=False)
     |
     |  Class summarizing bipartitions, branch lengths, and topologies from many trees
     |
     |  Method resolution order:
     |      BigTreeSummary
     |      TreeSummary
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(self, include_zeroterms=False, outgroup=None, rootmid=False)
     |      TreeSummary constructor. Initializes relevant data structures
     |
     |  add_tree(self, curtree, weight=1.0)
     |      Add tree to treesummary, update all summaries
     |
     |  topo_report(self)
     |      Returns list of [freq, treestring] lists
     |
     |  update(self, treesummary)
     |      Merge this class with other treesummary
     |
     |  ----------------------------------------------------------------------
     |  Methods inherited from TreeSummary:
     |
     |  bipart_report(self, includeleaves=True, minfreq=0.05)
     |      Return processed, almost directly printable, summary of all observed bipartitions
     |
     |  bipart_result(self)
     |      Return raw summary of all observed bipartitions
     |
     |  bipart_to_string(self, bipartition, position_dict)
     |      Takes bipartition (set of two leaf sets) and returns string representation
     |
     |  contree(self, cutoff=0.5, allcompat=False, lab='freq')
     |      Returns a consensus tree built from selected bipartitions
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from TreeSummary:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Branchstruct(builtins.object)
     |  Branchstruct(length=0.0, label='')
     |
     |  Class that emulates a struct. Keeps branch-related info
     |
     |  Methods defined here:
     |
     |  __init__(self, length=0.0, label='')
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Distmatrix(builtins.object)
     |  Class representing distance matrix for set of taxa. Knows hot to compute trees
     |
     |  Methods defined here:
     |
     |  __init__(self)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  __str__(self)
     |      Returns distance matrix as string
     |
     |  getdist(self, name1, name2)
     |      Returns distance between named entries
     |
     |  nj(self)
     |      Computes neighbor joining tree, returns Tree object
     |
     |  setdist(self, name1, name2, dist)
     |      Sets distance between named entries
     |
     |  upgma(self)
     |      Computes UPGMA tree, returns Tree object
     |
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |
     |  from_alignment(alignment, dist='pdist') from builtins.type
     |      Construct Distmatrix object from alignment object
     |
     |  from_distdict(distdict) from builtins.type
     |      Construct Distmatrix object from dictionary of {(name1, name2):dist}
     |
     |  from_numpy_array(nparray, namelist, n) from builtins.type
     |      Construct Distmatrix object from numpy array and corresponding list of names
     |
     |  from_phylip_string(dmat_string, n=None) from builtins.type
     |      Constructs Distmatrix object from string corresponding to PHYLIP square distance matrix
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Globals(builtins.object)
     |  Class containing globally used functions and labels.
     |
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  biparts = {}

    class Newicktreefile(Treefile)
     |  Newicktreefile(filename=None, filecontent=None)
     |
     |  Class representing Newick tree file. Iteration returns tree-objects
     |
     |  Method resolution order:
     |      Newicktreefile
     |      Treefile
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(self, filename=None, filecontent=None)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  __iter__(self)
     |
     |  __next__(self)
     |
     |  ----------------------------------------------------------------------
     |  Methods inherited from Treefile:
     |
     |  get_treestring(self)
     |      Return next tree-string
     |
     |  read_tree(self)
     |      Reads one tree from file and returns as Tree object. Returns None when exhausted file
     |
     |  read_trees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Treefile:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Nexustreefile(Treefile)
     |  Nexustreefile(filename=None, filecontent=None)
     |
     |  Class representing Nexus tree file. Iteration returns tree object or None
     |
     |  Method resolution order:
     |      Nexustreefile
     |      Treefile
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __init__(self, filename=None, filecontent=None)
     |      Read past NEXUS file header, parse translate block if present
     |
     |  __iter__(self)
     |
     |  __next__(self, noreturn=False)
     |
     |  ----------------------------------------------------------------------
     |  Methods inherited from Treefile:
     |
     |  get_treestring(self)
     |      Return next tree-string
     |
     |  read_tree(self)
     |      Reads one tree from file and returns as Tree object. Returns None when exhausted file
     |
     |  read_trees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Treefile:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Topostruct(builtins.object)
     |  Topostruct(count=1, treestring='')
     |
     |  Class that emulates a struct. Keeps topology-related info
     |
     |  Methods defined here:
     |
     |  __init__(self, count=1, treestring='')
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Tree(builtins.object)
     |  Class representing basic phylogenetic tree object.
     |
     |  Methods defined here:
     |
     |  __eq__(self, other, blenprecision=0.005)
     |      Implements equality testing for Tree objects
     |
     |  __hash__(self)
     |      Implements hashing for Tree objects, so they can be used as keys in dicts
     |
     |  __init__(self)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  __iter__(self)
     |      Returns iterator object for Tree object. Yields subtrees with .basalbranch attribute
     |
     |  __str__(self)
     |      Prints table of parent-child relationships including branch lengths and labels
     |
     |  add_branch(self, bipart, blen=0.0, label='')
     |      Adds branch represented by bipartition to unresolved tree.
     |
     |  average_ancdist(self, leaflist, return_median=False)
     |      Return average or median patristic distance from leaves to their MRCA
     |
     |  average_pairdist(self, leaflist, return_median=False)
     |      Return average or median pairwise, patristic distance between leaves in leaflist
     |
     |  bipdict(self)
     |      Returns tree in the form of a "bipartition dictionary"
     |
     |  build_dist_dict(self)
     |      Construct dictionary keeping track of all pairwise distances between nodes
     |
     |  build_parent_dict(self)
     |      Forces construction of parent_dict enabling faster lookups
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
     |      Divides tree into clusters by cutting across tree "cutoff" distance from root.
     |      Returns list containing sets of leafnames
     |
     |  cluster_n(self, nclust)
     |      Divides tree into 'nclust' clusters based on distance from root.
     |
     |      Returns tuple containing: list with sets of leafnames (one set per cluster)
     |                                list of basenodes of clusters
     |                                set containing all leaves that are put in clusters
     |
     |  collapse_clade(self, leaflist, newname='clade')
     |      Replaces clade (leaves in leaflist) with single leaf.
     |      Branch length is set to average dist from basenode parent to leaves
     |
     |  deroot(self)
     |      If root is at bifurcation: remove root node, connect adjacent nodes
     |
     |  diameter(self, return_leaves=False)
     |      Return diameter: longest leaf-leaf distance along tree.
     |      If return_leaves is True: Return tuple with (maxdist, Leaf1, Leaf2)
     |
     |  figtree(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6, colorlist=None, color='0000FF')
     |      Returns figtree format tree as a string
     |
     |  find_central_leaf(self, leaflist)
     |      Finds central leaf for the provided list of leaves.
     |      Defined as having approximately equal distance to the two farthest leaves in leaflist
     |
     |  find_common_leaf(self, leaflist)
     |      Finds common leaf for the provided list of leaves.
     |      Defined as having the smallest average distance to remaining leaves
     |
     |  find_most_distant(self, node1, nodeset)
     |      Finds node in nodeset that is most distant from node1
     |
     |  find_mrca(self, leafset)
     |      Finds Most Recent Common Ancestor for the provided set of leaves
     |
     |  findbasenode(self, leafset)
     |      Finds node that is at the base of all leaves in leafset.
     |
     |  getlabel(self, node1, node2)
     |      Gets label on branch connecting node1 and node2
     |
     |  graft(self, other, node1, node2=None, blen1=0, blen2=0, graftlabel=None)
     |      Graft other tree to self
     |
     |      tree2 (other) intnodes will be renamed if names clash with those in tree1.
     |      node1: node in tree1 (self) below which tree2 (other) will be grafted. Cannot be root1
     |      node2: node in tree2 (other) below which tree2 will be attached (default is root of tree2)
     |      blen1: length of branch added to tree1 below graftpoint (lower of two newly created branches)
     |      blen2: length of branch above graft point and below tree2 (upper of two newly created branches)
     |      graftlabel: prepend value of "label" to leaf names on t2 (e.g: "graft_s1")
     |
     |  height(self)
     |      Returns height of tree: Largest root-to-tip distance
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
     |      Returns list of leaf names sorted alphabetically
     |
     |  length(self)
     |      Returns tree length (sum of all branch lengths)
     |
     |  nameprune(self, sep='_', keep_pattern=None)
     |      Prune leaves based on name redundancy:
     |      Find subtrees where all leaves have same start of name (up to first "_")
     |
     |  nearest_n_leaves(self, leaf1, n_neighbors)
     |      Returns set of N leaves closest to leaf along tree (patristic distance)
     |
     |  nearleafs(self, leaf1, maxdist)
     |      Returns set of leaves that are less than maxdist from leaf, measured along branches
     |
     |  newick(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6)
     |      Returns Newick format tree string representation of tree object
     |
     |  nexus(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6)
     |      Returns nexus format tree as a string
     |
     |  nodedepth(self, node)
     |      Returns depth of node: distance from furthest leaf-level to node
     |
     |  nodedist(self, node1, node2=None)
     |      Returns distance between node1 and node2 along tree (patristic distance)
     |
     |  nodedistlist(self, node1, nodelist)
     |      Returns list of distances from node1 to nodes in nodelist (same order as nodelist)
     |
     |  nodepath(self, node1, node2)
     |      Returns path between node1 and node2 along tree.
     |
     |  nodepath_fromdict(self, node1, node2)
     |      Returns path between node1 and node2 along tree, from preconstructed path_dict
     |
     |  numberprune(self, nkeep, keeplist=None, keep_common_leaves=False, keep_most_distant=False, return_leaves=False, enforce_n=False)
     |      Prune tree so 'nkeep' leaves remain, approximately evenly spaced over tree.
     |
     |      "keeplist" can be used to specify leaves that _must_ be retained.
     |      'keep_common_leaves' requests preferential retainment of leaves with many neighbors
     |      (default is to keep leaves that are as equally spaced as possible)
     |      'keep_most_distant' requests that the two most distant leaves in tree
     |                  (which spread out the diameter) should be kept
     |      'return_leaves': return selected leaves, but do not actually prune tree
     |      'enforce_n' enforce exactly N leaves in pruned tree
     |                  (normally leaves in includelist and most distant are additional to N)
     |
     |  parent(self, node)
     |      Returns parent of node
     |
     |  prune_maxlen(self, nkeep, return_leaves=False)
     |      Prune tree so remaining nkeep leaves spread out maximal percentage of branch length
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
     |      Removes leaves in list from tree, cleans up so remaining tree structure is sane
     |
     |  rename_intnode(self, oldnum, newnum)
     |      Changes number of one internal node
     |
     |  rename_leaf(self, oldname, newname, fixdups=False)
     |      Changes name of one leaf. Automatically fixes duplicates if requested
     |
     |  reroot(self, node1, node2=None, polytomy=False, node1dist=0.0)
     |      Places new root on branch between node1 and node2, node1dist from node1
     |
     |  resolve(self)
     |      Randomly resolves multifurcating tree by by adding zero-length internal branches.
     |
     |  rootmid(self)
     |      Performs midpoint rooting of tree
     |
     |  rootminvar(self)
     |      Performs minimum variance rooting of tree
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
     |      Returns sorted intnode list for breadth-first traversal of tree
     |
     |  spr(self, subtree_node, regraft_node)
     |      Subtree Pruning and Regrafting.
     |
     |      subtree_node: basenode of subtree that will be pruned.
     |      regraft_node: node in tree below which subtree will be grafted
     |
     |  subtree(self, basenode, return_basalbranch=False)
     |      Returns subtree rooted at basenode as Tree object
     |
     |  topology(self)
     |      Returns set of sets of sets representation of topology ("naked bipdict")
     |
     |  transname(self, namefile)
     |      Translate all leaf names using oldname/newname pairs in namefile
     |
     |  treedist(self, other, normalise=True, verbose=False)
     |      Compute symmetric tree distance (Robinson Foulds) between self and other tree.
     |      Normalised measure returned by default
     |
     |  treesim(self, other, verbose=False)
     |      Compute normalised symmetric similarity between self and other tree
     |
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |
     |  from_biplist(biplist) from builtins.type
     |      Constructor 2: Tree object from bipartition list
     |
     |  from_leaves(leaflist) from builtins.type
     |      Constructor 3: star-tree object from list of leaves
     |
     |  from_string(orig_treestring, transdict=None) from builtins.type
     |      Constructor 1: Tree object from tree-string in Newick format
     |
     |  randtree(leaflist=None, ntips=None, randomlen=False, name_prefix='s') from builtins.type
     |      Constructor 4: tree with random topology from list of leaf names OR number of tips
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class TreeError(builtins.Exception)
     |  Method resolution order:
     |      TreeError
     |      builtins.Exception
     |      builtins.BaseException
     |      builtins.object
     |
     |  Data descriptors defined here:
     |
     |  __weakref__
     |      list of weak references to the object (if defined)
     |
     |  ----------------------------------------------------------------------
     |  Methods inherited from builtins.Exception:
     |
     |  __init__(self, /, *args, **kwargs)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  ----------------------------------------------------------------------
     |  Static methods inherited from builtins.Exception:
     |
     |  __new__(*args, **kwargs) from builtins.type
     |      Create and return a new object.  See help(type) for accurate signature.
     |
     |  ----------------------------------------------------------------------
     |  Methods inherited from builtins.BaseException:
     |
     |  __delattr__(self, name, /)
     |      Implement delattr(self, name).
     |
     |  __getattribute__(self, name, /)
     |      Return getattr(self, name).
     |
     |  __reduce__(...)
     |      Helper for pickle.
     |
     |  __repr__(self, /)
     |      Return repr(self).
     |
     |  __setattr__(self, name, value, /)
     |      Implement setattr(self, name, value).
     |
     |  __setstate__(...)
     |
     |  __str__(self, /)
     |      Return str(self).
     |
     |  with_traceback(...)
     |      Exception.with_traceback(tb) --
     |      set self.__traceback__ to tb and return self.
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from builtins.BaseException:
     |
     |  __cause__
     |      exception cause
     |
     |  __context__
     |      exception context
     |
     |  __dict__
     |
     |  __suppress_context__
     |
     |  __traceback__
     |
     |  args

    class TreeSet(builtins.object)
     |  Class for storing and manipulating a number of trees
     |
     |  Methods defined here:
     |
     |  __getitem__(self, index)
     |      Implements indexing of treeset.
     |
     |      Simple index returns single tree.
     |      Slice returns TreeSet object with selected subset of trees
     |
     |  __init__(self)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  __iter__(self)
     |      Returns fresh iterator object allowing iteration over Treeset (which is itself an iterable)
     |
     |  __len__(self)
     |
     |  addtree(self, tree)
     |      Adds Tree object to Treeset object
     |
     |  addtreeset(self, treeset)
     |      Adds all trees in TreeSet object to this TreeSet object
     |
     |  newick(self, printdist=True, printlabels=True)
     |      Returns newick format tree as a string
     |
     |  nexus(self, printdist=True, printlabels=True)
     |      Returns nexus format tree as a string
     |
     |  rootmid(self)
     |      Performs midpoint rooting on all trees in TreeSet
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  TreeSetIterator = <class 'phylotreelib.TreeSet.TreeSetIterator'>

    class TreeSummary(builtins.object)
     |  TreeSummary(include_zeroterms=False)
     |
     |  Class summarizing bipartitions and branch lengths (but not topologies) from many trees
     |
     |  Methods defined here:
     |
     |  __init__(self, include_zeroterms=False)
     |      TreeSummary constructor. Initializes relevant data structures
     |
     |  add_tree(self, curtree, weight=1.0)
     |      Add tree object to treesummary, update all relevant summaries
     |
     |  bipart_report(self, includeleaves=True, minfreq=0.05)
     |      Return processed, almost directly printable, summary of all observed bipartitions
     |
     |  bipart_result(self)
     |      Return raw summary of all observed bipartitions
     |
     |  bipart_to_string(self, bipartition, position_dict)
     |      Takes bipartition (set of two leaf sets) and returns string representation
     |
     |  contree(self, cutoff=0.5, allcompat=False, lab='freq')
     |      Returns a consensus tree built from selected bipartitions
     |
     |  update(self, treesummary)
     |      Merge this class with external treesummary
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Treefile(builtins.object)
     |  Treefile(filename=None, filecontent=None)
     |
     |  Abstract base-class for representing tree file objects.
     |
     |  Methods defined here:
     |
     |  __init__(self, filename=None, filecontent=None)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  get_treestring(self)
     |      Return next tree-string
     |
     |  read_tree(self)
     |      Reads one tree from file and returns as Tree object. Returns None when exhausted file
     |
     |  read_trees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

FUNCTIONS
    escape_metachars(text, metachars='.^$*+?{}[]\\|()')
        Add backslashes to escape metachars in input string.

    main()
        # # Placeholder: Insert test code here and run module in standalone mode

    remove_comments(text, leftdelim, rightdelim=None)
        Takes input string and strips away commented text, delimited by 'leftdelim' and 'rightdelim'.
        Also deals with nested comments.
```


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
        Interner
        Topostruct
        Tree
        TreeSet
        TreeSummary
            BigTreeSummary
        Treefile
        TreefileBase
            Newicktreefile
            Nexustreefile
    
    class BigTreeSummary(TreeSummary)
     |  BigTreeSummary(interner=None, store_trees=False)
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
     |  __init__(self, interner=None, store_trees=False)
     |      TreeSummary constructor. Initializes relevant data structures
     |  
     |  add_tree(self, curtree, weight=1.0)
     |      Add tree to treesummary, update all summaries
     |  
     |  max_clade_cred_tree(self, labeldigits=3)
     |      Find and return maximum clade credibility tree
     |  
     |  update(self, other)
     |      Merge this object with other treesummary
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  toposummary
     |      Property method for lazy evaluation of topostruct.freq
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from TreeSummary:
     |  
     |  __len__(self)
     |  
     |  add_branchid(self)
     |      Adds attribute .branchID to all bipartitions in .bipartsummary
     |      External bipartitions are labeled with the leafname.
     |      Internal bipartitions are labeled with consecutive numbers by decreasing frequency
     |  
     |  contree(self, cutoff=0.5, allcompat=False, labeldigits=3)
     |      Returns a consensus tree built from selected bipartitions
     |  
     |  log_clade_credibility(self, topology)
     |      Compute log clade credibility for topology (sum of log(freq) for all branches)
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties inherited from TreeSummary:
     |  
     |  bipartsummary
     |      Property method for lazy evaluation of freq, var, and sem for branches
     |  
     |  sorted_biplist
     |      Return list of bipartitions.
     |      First external (leaf) bipartitions sorted by leafname.
     |      Then internal bipartitions sorted by freq
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
     |  copy(self)
     |      Returns copy of Branchstruct object, with all attributes included
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
     |  Class representing distance matrix for set of taxa. Knows how to compute trees
     |  
     |  Methods defined here:
     |  
     |  __init__(self)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  __str__(self)
     |      Returns distance matrix as string
     |  
     |  avdist(self)
     |      Returns average dist in matrix (not including diagonal)
     |  
     |  clean_names(self, illegal=',:;()[]', rep='_')
     |      Rename items to avoid characters that are problematic in Newick tree strings:
     |      Replaces all occurrences of chars in 'illegal' by 'rep'
     |  
     |  getdist(self, name1, name2)
     |      Returns distance between named entries
     |  
     |  nj(self)
     |      Computes neighbor joining tree, returns Tree object
     |  
     |  rename(self, oldname, newname)
     |      Changes name of one item from oldname to newname
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
     |  from_distdict(distdict) from builtins.type
     |      Construct Distmatrix object from nested dictionary of dists: distdict[name1][name2] = dist
     |  
     |  from_distfile(distfilename) from builtins.type
     |      Construct Distmatrix object from file containing rows of: name1 name2 distance
     |  
     |  from_numpy_array(nparray, namelist) from builtins.type
     |      Construct Distmatrix object from numpy array and corresponding list of names
     |      Names in namelist must be in same order as indices in numpy 2D array
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
    
    class Interner(builtins.object)
     |  Class used for interning various objects.
     |  
     |  Methods defined here:
     |  
     |  __init__(self)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  intern_bipart(self, bipart)
     |  
     |  intern_leafset(self, leafset)
     |  
     |  intern_topology(self, topology)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class Newicktreefile(TreefileBase)
     |  Newicktreefile(filename=None, filecontent=None)
     |  
     |  Class representing Newick tree file. Iteration returns tree-objects
     |  
     |  Method resolution order:
     |      Newicktreefile
     |      TreefileBase
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
     |  Methods inherited from TreefileBase:
     |  
     |  __enter__(self)
     |      Implements context manager behaviour for TreefileBase types.
     |      Usage example:
     |          with pt.Newicktreefile(filename) as tf:
     |              mytree = tf.readtree()
     |          mytree.rootminvar()
     |  
     |  __exit__(self, type, value, traceback)
     |      Implements context manager behaviour for TreefileBase types.
     |      Usage example:
     |          with pt.Newicktreefile(filename) as tf:
     |              mytree = tf.readtree()
     |          mytree.rootminvar()
     |  
     |  close(self)
     |      For explicit closing of Treefile before content exhausted
     |  
     |  get_treestring(self)
     |      Return next tree-string
     |  
     |  readtree(self)
     |      Reads one tree from file and returns as Tree object. Returns None when exhausted file
     |  
     |  readtrees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from TreefileBase:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class Nexustreefile(TreefileBase)
     |  Nexustreefile(filename=None, filecontent=None)
     |  
     |  Class representing Nexus tree file. Iteration returns tree object or None
     |  
     |  Method resolution order:
     |      Nexustreefile
     |      TreefileBase
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
     |  Methods inherited from TreefileBase:
     |  
     |  __enter__(self)
     |      Implements context manager behaviour for TreefileBase types.
     |      Usage example:
     |          with pt.Newicktreefile(filename) as tf:
     |              mytree = tf.readtree()
     |          mytree.rootminvar()
     |  
     |  __exit__(self, type, value, traceback)
     |      Implements context manager behaviour for TreefileBase types.
     |      Usage example:
     |          with pt.Newicktreefile(filename) as tf:
     |              mytree = tf.readtree()
     |          mytree.rootminvar()
     |  
     |  close(self)
     |      For explicit closing of Treefile before content exhausted
     |  
     |  get_treestring(self)
     |      Return next tree-string
     |  
     |  readtree(self)
     |      Reads one tree from file and returns as Tree object. Returns None when exhausted file
     |  
     |  readtrees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from TreefileBase:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class Topostruct(builtins.object)
     |  Class that emulates a struct. Keeps topology-related info
     |  
     |  Data descriptors defined here:
     |  
     |  freq
     |  
     |  tree
     |  
     |  weight
    
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
     |  add_branch(self, bipart, branchstruct)
     |      Adds branch represented by bipartition to unresolved tree.
     |  
     |  add_leaf(self, parent, newleafname, branchstruct)
     |      Adds new leaf to existing intnode ´parent´
     |  
     |  average_ancdist(self, leaflist, return_median=False)
     |      Return average or median patristic distance from leaves to their MRCA
     |  
     |  average_pairdist(self, leaflist, return_median=False)
     |      Return average or median pairwise, patristic distance between leaves in leaflist
     |  
     |  bipdict(self, interner=None)
     |      Returns tree in the form of a "bipartition dictionary"
     |  
     |  build_dist_dict(self)
     |      Construct dictionary keeping track of all pairwise distances between nodes
     |  
     |  build_parent_dict(self)
     |      Constructs _parent_dict enabling faster lookups, when needed
     |  
     |  build_path_dict(self)
     |      Construct dictionary keeping track of all pairwise paths between nodes
     |  
     |  check_bip_compatibility(self, bipart)
     |      Checks the compatibility between bipartition and tree.
     |      Returns tuple of: is_present, is_compatible, insert_tuple
     |              where insert_tuple = None or (parentnode, childmovelist)
     |      is_present:
     |          True if bipartition is already present in tree. Implies "is_compatible = True"
     |      is_compatible:
     |          True if bipartition is compatible with tree. "is_present" can be True or False
     |      insert_tuple:
     |          If is_compatible: Tuple of (parentnode, childmovelist) parameters for insert_node
     |          If not is_compatible: None
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
     |  
     |  collapse_clade(self, leaflist, newname='clade')
     |      Replaces clade (leaves in leaflist) with single leaf.
     |      Branch length is set to average dist from basenode parent to leaves
     |  
     |  copy_treeobject(self, copylengths=True, copylabels=True)
     |      Returns copy of Tree object. Copies structure and branch lengths.
     |      Caches and any user-added attributes are not copied.
     |      Similar to effect of copy.deepcopy but customized and much faster
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
     |  find_mrca(self, leaves)
     |      Finds Most Recent Common Ancestor for the provided set of leaves
     |  
     |  findbasenode(self, leafset)
     |      Finds node that is at the base of all leaves in leafset.
     |  
     |  get_branchstruct(self, node1, node2)
     |      Returns Branchstruct object from branch between node1 and node2
     |  
     |  getlabel(self, node1, node2)
     |      Gets label on branch connecting node1 and node2
     |  
     |  graft(self, other, node1, node2=None, blen1=0, blen2=0, graftlabel=None, graft_with_other_root=False)
     |      Graft other tree to self
     |      
     |      tree2 (other) intnodes will be renamed if names clash with those in tree1.
     |      node1: node in tree1 (self) below which tree2 (other) will be grafted. Cannot be root1
     |      node2: node in tree2 (other) below which tree2 will be attached (default is root of tree2)
     |      blen1: length of branch added to tree1 below graftpoint (lower of two newly created branches)
     |      blen2: length of branch above graft point and below tree2 (upper of two newly created branches)
     |      graftlabel: prepend value of "label" to leaf names on t2 (e.g: "graft_s1")
     |      graft_with_other_root: use root of other as graftpoint (i.e., do not add extra basal
     |                             branch between other.root and self.graftpoint)
     |  
     |  has_same_root(self, other)
     |      Compares two trees. Returns True if topologies are same and rooted in same place
     |  
     |  height(self)
     |      Returns height of tree: Largest root-to-tip distance
     |  
     |  insert_node(self, parent, childnodes, branchstruct)
     |      Inserts an extra node between parent and children listed in childnodes list
     |      (so childnodes are now attached to newnode instead of parent).
     |      The branchstruct will be attached to the branch between parent and newnode.
     |      Branches to childnodes retain their original branchstructs.
     |      The node number of the new node is returned
     |  
     |  is_bifurcation(self, node)
     |      Checks if internal node is at bifurcation (has two children)
     |  
     |  is_compatible_with(self, bipart)
     |      Checks whether a given bipartition is compatible with the tree.
     |      Note: also returns True if bipartition is already in tree
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
     |  match_nodes(self, other)
     |      Compares two identical trees with potentially different internal node IDs.
     |      Returns tuple containing following:
     |          Dictionary giving mapping from nodeid in self to nodeid in other (also leaves)
     |          unmatched_root1: "None" or id of unmatched root in self if root at bifurcation
     |          unmatched_root2: "None" or id of unmatched root in other if root at bifurcation
     |      
     |      Note: The last two are only different from None if the trees dont have the same
     |      exact rooting
     |  
     |  n_bipartitions(self)
     |      Returns the number of bipartitions (= number of internal branches) in tree
     |      Note: if root is at bifurcation, then those 2 branches = 1 bipartition
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
     |  newick(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6, labelfield='label', transdict=None)
     |      Returns Newick format tree string representation of tree object
     |  
     |  nexus(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6, labelfield='label', translateblock=False)
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
     |  patristic_distdict(self)
     |      Return nested dictionary giving all pairwise, patristic distances:
     |      dict[node1][node2] = patristic distance
     |  
     |  prune_maxlen(self, nkeep, return_leaves=False)
     |      Prune tree so remaining nkeep leaves spread out maximal percentage of branch length
     |  
     |  remote_children(self, parent)
     |      Returns set containing all leaves that are descendants of parent
     |  
     |  remote_nodes(self, parent)
     |      Returns set containing all nodes (intnodes and leaves) that are descendants of parent.
     |      This set includes parent itself
     |  
     |  remove_branch(self, node1, node2)
     |      Removes branch connecting node1 and node2 (thereby possibly creating polytomy)
     |      Length of removed branch is distributed among descendant branches.
     |      This means tree length is conserved.
     |      Descendant nodes will be farther apart from each other, but closer to outside nodes.
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
     |  set_branch_attribute(self, node1, node2, attrname, attrvalue)
     |      Set the value of any branch attribute.
     |      attrname: Name of attribute (e.g., "length")
     |      attrvalue: Value of attribute (e.g. 0.153)
     |  
     |  set_nodeid_labels(self)
     |      Sets labels to be the same as the child node ID
     |      Allows use of e.g. Figtree to show nodeIDs as nodelabels
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
     |  spr(self, subtree_node=None, regraft_node=None)
     |      Subtree Pruning and Regrafting.
     |      
     |      subtree_node: basenode of subtree that will be pruned.
     |      regraft_node: node in tree below which subtree will be grafted
     |      
     |      If no parameters are specified (both are None): perform random SPR
     |      If only subtree_node is specified: choose random regraft_node
     |      
     |      Must specify either both parameters, no parameters, or only subtree_node
     |  
     |  subtree(self, basenode, return_basalbranch=False)
     |      Returns subtree rooted at basenode as Tree object
     |  
     |  topology(self)
     |      Returns set of sets of sets representation of topology ("naked bipdict")
     |  
     |  transdict(self)
     |      Returns dictionary of {name:number_as_string} for use in translateblocks
     |  
     |  translateblock(self, transdict)
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
     |      Constructor: Tree object from bipartition list
     |  
     |  from_branchinfo(parentlist, childlist, lenlist=None, lablist=None) from builtins.type
     |      Constructor: Tree object from information about all branches in tree
     |      
     |      Information about one branch is conceptually given as:
     |          parentnodeID, childnodeID, [length], [label]
     |      
     |      The function takes as input 2 to 4 separate lists containing:
     |          IDs of parents (internal nodes, so integer values)
     |          ID of children (internal or leaf nodes, so integer or string)
     |          Length of branches (optional)
     |          Label of branches (optional)
     |      
     |      The four lists are assumed to have same length and be in same order (so index n in
     |      each list corresponds to same branch).
     |      
     |      Note: most IDs appear multiple times in lists
     |      Note 2: can be used as workaround so user can specify IDs for internal nodes
     |  
     |  from_leaves(leaflist) from builtins.type
     |      Constructor: star-tree object from list of leaves
     |  
     |  from_string(orig_treestring, transdict=None) from builtins.type
     |      Constructor: Tree object from tree-string in Newick format
     |  
     |  from_topology(topology) from builtins.type
     |      Constructor: Tree object from topology
     |  
     |  randtree(leaflist=None, ntips=None, randomlen=False, name_prefix='s') from builtins.type
     |      Constructor: tree with random topology from list of leaf names OR number of tips
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  parent_dict
     |      Lazy evaluation of _parent_dict when needed
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
     |  Class for storing and manipulating a number of trees, which all have the same leafs
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
     |  nexus(self, printdist=True, printlabels=True, translateblock=True)
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
     |  TreeSummary(interner=None)
     |  
     |  Class summarizing bipartitions and branch lengths (but not topologies) from many trees
     |  
     |  Methods defined here:
     |  
     |  __init__(self, interner=None)
     |      TreeSummary constructor. Initializes relevant data structures
     |  
     |  __len__(self)
     |  
     |  add_branchid(self)
     |      Adds attribute .branchID to all bipartitions in .bipartsummary
     |      External bipartitions are labeled with the leafname.
     |      Internal bipartitions are labeled with consecutive numbers by decreasing frequency
     |  
     |  add_tree(self, curtree, weight=1.0)
     |      Add tree object to treesummary, update all relevant bipartition summaries
     |  
     |  contree(self, cutoff=0.5, allcompat=False, labeldigits=3)
     |      Returns a consensus tree built from selected bipartitions
     |  
     |  log_clade_credibility(self, topology)
     |      Compute log clade credibility for topology (sum of log(freq) for all branches)
     |  
     |  max_clade_cred_tree(self, filelist, skiplist=None, labeldigits=3)
     |      Find and return maximum clade credibility tree.
     |      Note: this version based on external treefile (and bipartsummary).
     |      Skiplist possibly contains number of trees to skip in each file (burnin)
     |  
     |  update(self, other)
     |      Merge this object with external treesummary
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  bipartsummary
     |      Property method for lazy evaluation of freq, var, and sem for branches
     |  
     |  sorted_biplist
     |      Return list of bipartitions.
     |      First external (leaf) bipartitions sorted by leafname.
     |      Then internal bipartitions sorted by freq
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
     |  Treefile(filename)
     |  
     |  Factory for making Newick or Nexus treefile objects. Autodetects fileformat
     |  
     |  Static methods defined here:
     |  
     |  __new__(klass, filename)
     |      Create and return a new object.  See help(type) for accurate signature.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class TreefileBase(builtins.object)
     |  TreefileBase(filename=None, filecontent=None)
     |  
     |  Abstract base-class for representing tree file objects.
     |  
     |  Methods defined here:
     |  
     |  __enter__(self)
     |      Implements context manager behaviour for TreefileBase types.
     |      Usage example:
     |          with pt.Newicktreefile(filename) as tf:
     |              mytree = tf.readtree()
     |          mytree.rootminvar()
     |  
     |  __exit__(self, type, value, traceback)
     |      Implements context manager behaviour for TreefileBase types.
     |      Usage example:
     |          with pt.Newicktreefile(filename) as tf:
     |              mytree = tf.readtree()
     |          mytree.rootminvar()
     |  
     |  __init__(self, filename=None, filecontent=None)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  close(self)
     |      For explicit closing of Treefile before content exhausted
     |  
     |  get_treestring(self)
     |      Return next tree-string
     |  
     |  readtree(self)
     |      Reads one tree from file and returns as Tree object. Returns None when exhausted file
     |  
     |  readtrees(self, discardprop=0.0)
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
    main()
        # # Placeholder: Insert test code here and run module in standalone mode
    
    remove_comments(text, leftdelim, rightdelim=None)
        Takes input string and strips away commented text, delimited by 'leftdelim' and 'rightdelim'.
        Also deals with nested comments.

FILE
    /Users/agpe/Documents/code/python/phylotreelib/phylotreelib.py



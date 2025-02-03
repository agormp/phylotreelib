Help on module phylotreelib:

NAME
    phylotreelib - Classes and methods for analyzing, manipulating, and building phylogenetic trees

CLASSES
    builtins.Exception(builtins.BaseException)
        TreeError
    builtins.object
        Bipartition
        Branchstruct
        Clade
        Distmatrix
        Interner
        NewickStringParser
        Nodestruct
        RootBipStruct
        SymmetricMatrix
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
     |  BigTreeSummary(store_trees=False, trackbips=True, trackclades=False, trackroot=False, trackblen=False, trackdepth=False)
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
     |  __init__(self, store_trees=False, trackbips=True, trackclades=False, trackroot=False, trackblen=False, trackdepth=False)
     |      TreeSummary constructor. Initializes relevant data structures
     |
     |  add_tree(self, curtree, weight=1.0)
     |      Add tree to treesummary, update all summaries
     |
     |  max_bipart_cred_tree(self)
     |      Find maximum bipartition credibility tree. Return tuple of (maxcredtree, maxlogcred)
     |
     |  max_clade_cred_tree(self)
     |      Find maximum clade credibility tree. Return tuple of (maxcredtree, maxlogcred)
     |
     |  update(self, other)
     |      Merge this object with other treesummary
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  biptoposummary
     |      Property method for lazy evaluation of topostruct.posterior
     |
     |  cladetoposummary
     |      Property method for lazy evaluation of topostruct.posterior
     |
     |  ----------------------------------------------------------------------
     |  Methods inherited from TreeSummary:
     |
     |  __len__(self)
     |
     |  compute_rootcred(self, tree)
     |      Returns root credibility (frequency of tree's root among observed trees) based
     |      on current root of sum_tree and information in self._rootbip_summary
     |
     |  contree(self, cutoff=0.5, allcompat=False)
     |      Returns a consensus tree built from selected bipartitions.
     |
     |  log_bipart_credibility(self, biptopology)
     |      Compute log bipartition credibility for topology (sum of log(freq) for all branches)
     |
     |  log_clade_credibility(self, cladetopology)
     |      Compute log clade credibility for topology (sum of log(freq) for all clades)
     |
     |  root_maxfreq(self, sum_tree)
     |      Uses info about root bipartitions in TreeSummary to place root on summary tree.
     |      Also sets tree attribute rootcred:
     |          probability (freq among input trees) of current location of root
     |
     |      If tree has branch lengths:
     |      Divides length of root bipartition among two branches in accordance with average
     |      fraction of lengths seen for this rootbip across all trees.
     |
     |  set_ca_node_depths(self, sum_tree, wt_count_burnin_filename_list)
     |      Set branch lengths on summary tree based on mean node depth for clades corresponding
     |      to MRCA of clade's leaves. (same as "--height ca" in BEAST's treeannotator)
     |      This means that all input trees are used when computing
     |      mean for each node (not just the input trees where that exact monophyletic clade
     |      is present)
     |
     |  set_clade_credibility(self, tree, precision=6)
     |      Set clade credibility on provided target tree based on freq of clade in TreeSummary.
     |
     |      NOTE: only works if all clades in tree have been observed at least once. The option
     |              will therefore not work with all rootings
     |
     |  set_mean_biplen(self, sum_tree)
     |      Sets branch-length, -var, and -sem for each branch in sum_tree,
     |      based on values in bipsummary.
     |      Should only be called when sum_tree constructed based on clades (MCC),
     |      but brlens are to be set based on mean bipartition values
     |
     |  set_mean_node_depths(self, sum_tree)
     |      Set branch lengths on summary tree based on mean node depth for clades corresponding
     |      to parent and child nodes (blen = depth_parent - depth_child).
     |
     |      NOTE 1: only meaningful if input trees are based on a clock model.
     |      NOTE 2: only works if all clades in tree have been observed at least once. The option
     |              will therefore not work with all rootings, and may also fail for majority rule
     |              consensus trees
     |      NOTE 3: only uses node depths from monophyletic clades (so some values may be set
     |      based on very few trees)
     |
     |  set_rootcredibility(self, sum_tree, precision=6)
     |      Returns sum_tree with root credibilities as attributes on each branch
     |      rootcred = fraction of trees in input set where the root was on this branch (bipartition)
     |      If root was never on a branch: assign the value 0.0
     |      Added as attribute .rootcred to Branchstruct for branches on sum_tree
     |      Also sets tree-attribute cumrootcred:
     |          sum of rootcredibilities for all branches (bipartitions) included on tree
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties inherited from TreeSummary:
     |
     |  bipartsummary
     |      Property method for lazy evaluation of freq, var, sd, and sem for branches
     |
     |  cladesummary
     |      Property method for lazy evaluation of freq, var, and sem for node depths
     |
     |  rootbipsummary
     |      Property method for lazy evaluation of freq (=rootcred) for rootbips
     |
     |  sorted_biplist
     |      Return list of bipartitions.
     |      First external (leaf) bipartitions sorted by leafname.
     |      Then internal bipartitions sorted by freq
     |
     |  sorted_rootbips
     |      Return list of root-bipartitions (branches where root has been seen), sorted by
     |      occurrence (count on tree samples added to TreeSummary)
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from TreeSummary:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class Bipartition(builtins.object)
     |  Bipartition(leafset1, all_leaves_set, sorted_leaf_list, leaf2index)
     |
     |  Class that represents a bipartition: the split of leaves in two sets that corresponds
     |  to a branch in a tree. Class has a fast method for hashing and therefore useful when
     |  comparing Bipartitions
     |
     |  Methods defined here:
     |
     |  __eq__(self, other)
     |      Return self==value.
     |
     |  __hash__(self)
     |      Return hash(self).
     |
     |  __init__(self, leafset1, all_leaves_set, sorted_leaf_list, leaf2index)
     |      Initialise bipartition objects based on only one half of bipartition.
     |      If the complement of leafset1 is smaller, then that will be stored instead.
     |      If they have same size: store leafset with smaller hash.
     |
     |      leafset1: one half of bipartition (to save time in caller)
     |      all_leaves_set: set of all leaf names
     |      sorted_leaf_list: sorted list of all leaf names
     |      leaf2index: dict of {leaf:index in sorted_leaf_list}
     |
     |      Last 3 params will typically point to attributes in parent Tree object
     |
     |  __iter__(self)
     |      # Python note: this allows unpacking as if the class was a tuple: bip1, bip2 = bipartition
     |
     |  __repr__(self)
     |      Return repr(self).
     |
     |  __str__(self)
     |      Return str(self).
     |
     |  get_bipartitions(self)
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  indices
     |
     |  leaf2index
     |
     |  leaf_list
     |
     |  leaf_set

    class Branchstruct(builtins.object)
     |  Branchstruct(length=0.0, **attributes)
     |
     |  Class that emulates a struct. Keeps branch-related info
     |
     |  Methods defined here:
     |
     |  __init__(self, length=0.0, **attributes)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  __repr__(self)
     |      Return repr(self).
     |
     |  __str__(self)
     |      Return str(self).
     |
     |  copy(self)
     |      Returns copy of Branchstruct object, with all attributes included
     |
     |  merge(self, other, check_compat=False)
     |      Merges two Branchstructs and returns new Branchstruct.
     |      Useful for collapsing root branches.
     |      Branch lengths are summed, other attributes are taken from self
     |      check_compat: if True, check that non-length attributes in self and other match,
     |                    raise TreeError if not
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class Clade(builtins.object)
     |  Clade(leafset, all_leaves_set, sorted_leaf_list, leaf2index)
     |
     |  Class that represents a clade: the set of leaves descended from an internal node
     |
     |  Methods defined here:
     |
     |  __eq__(self, other)
     |      Return self==value.
     |
     |  __hash__(self)
     |      Return hash(self).
     |
     |  __init__(self, leafset, all_leaves_set, sorted_leaf_list, leaf2index)
     |      Initialise clade objects.
     |
     |      leafset: set of leaves descending from an internal node
     |      all_leaves_set: set of all leaf names in tree
     |      sorted_leaf_list: sorted list of all leaf names in tree
     |      leaf2index: dict of {leaf:index in sorted_leaf_list}
     |
     |      Last 3 params will typically point to attributes in parent Tree object
     |
     |  __iter__(self)
     |      # Python note: this allows unpacking as if the class was a tuple: c1 = myclade
     |
     |  __repr__(self)
     |      Return repr(self).
     |
     |  __str__(self)
     |      Return str(self).
     |
     |  get_clade(self)
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  all_leaves_set
     |
     |  indices
     |
     |  leaf2index
     |
     |  leaf_list

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
     |  from_distdict(distdict)
     |      Construct Distmatrix object from nested dictionary of dists: distdict[name1][name2] = dist
     |
     |  from_distfile(distfilename)
     |      Construct Distmatrix object from file containing rows of: name1 name2 distance
     |
     |  from_numpy_array(nparray, namelist)
     |      Construct Distmatrix object from numpy array and corresponding list of names
     |      Names in namelist must be in same order as indices in numpy 2D array
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

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
     |  intern_clade(self, clade)
     |
     |  intern_leafset(self, leafset)
     |
     |  store_unhashable(self, name, obj)
     |      Stores the unhashable object only if it hasn't been stored before.
     |      Always returns the stored object (whether newly stored or previously stored).
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class NewickStringParser(builtins.object)
     |  NewickStringParser(transdict=None, label_attr_name='label', label_type=<class 'str'>)
     |
     |  Class creating parser for specific Newick tree string. Used in Tree.from_string
     |
     |  Methods defined here:
     |
     |  __init__(self, transdict=None, label_attr_name='label', label_type=<class 'str'>)
     |      Initialise parser object.
     |      label_attr_name: name of label attribute on Branchstructs
     |      label_type: type conversion to perform on input label string (e.g. float)
     |
     |      To create derived class:
     |        (1) Call super __init__ in __init__ of derived class
     |        (2) Override set_token_delimiters() (return string of single-char delimiters)
     |        (3) Override update_dispatch() to modify self.dispatch dict
     |
     |  create_compiled_regex(self, delimiters)
     |      Creates compiled regex for faster parsing
     |
     |  parse(self, treeobj, treestring)
     |
     |  sanitychecks(self, treeobj, treestring)
     |
     |  set_token_delimiters(self)
     |      Override to modify the set of delimiters used when parsing.
     |      All delimiters must be single-character (handle multi-char delims in parsing code)
     |      Specify single string of delimiters that is returned
     |
     |  update_dispatch(self)
     |      Override to modify dispatch dictionary (extra edges and nodes in state diagram)
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class Newicktreefile(TreefileBase)
     |  Newicktreefile(filename=None, filecontent=None, interner=None)
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
     |  __init__(self, filename=None, filecontent=None, interner=None)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  __iter__(self)
     |
     |  __next__(self, returntree=True)
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
     |  readtree(self, returntree=True)
     |      Reads one treestring from file and returns as Tree object if requested.
     |      Returns None when exhausted file
     |
     |  readtrees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from TreefileBase:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class Nexustreefile(TreefileBase)
     |  Nexustreefile(filename=None, filecontent=None, interner=None)
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
     |  __init__(self, filename=None, filecontent=None, interner=None)
     |      Read past NEXUS file header, parse translate block if present
     |
     |  __iter__(self)
     |
     |  __next__(self, returntree=True)
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
     |  readtree(self, returntree=True)
     |      Reads one treestring from file and returns as Tree object if requested.
     |      Returns None when exhausted file
     |
     |  readtrees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from TreefileBase:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class Nodestruct(builtins.object)
     |  Nodestruct(depth=0.0)
     |
     |  Class that emulates a struct. Keeps node-related info
     |
     |  Methods defined here:
     |
     |  __init__(self, depth=0.0)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  __repr__(self)
     |      Return repr(self).
     |
     |  __str__(self)
     |      Return str(self).
     |
     |  copy(self)
     |      Returns copy of Nodestruct object, with all attributes included
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class RootBipStruct(builtins.object)
     |  RootBipStruct(leafset1, blen1, leafset2, blen2)
     |
     |  Methods defined here:
     |
     |  __init__(self, leafset1, blen1, leafset2, blen2)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  add(self, leafset1, blen1, leafset2, blen2)
     |      Adds branch length fractions to the current sum and increments the count.
     |
     |  avg_frac(self, leafset)
     |
     |  merge(self, other)
     |      Merges this RootBipStruct with another (for same bipartition)
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class SymmetricMatrix(builtins.object)
     |  Base class for symmetric matrices (especially distance matrices for trees or taxa).
     |  Only needs to set upper half, but can be accessed using any order of indices
     |
     |  Methods defined here:
     |
     |  __init__(self)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class Topostruct(builtins.object)
     |  Class that emulates a struct. Keeps topology-related info
     |
     |  Data descriptors defined here:
     |
     |  posterior
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
     |      Prints table of parent-child relationships including branch lengths and other attributes.
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
     |  bipart_is_present(self, bipart)
     |      Checks whether a given bipartition is present in tree
     |
     |  bipdict(self, keep_remchild_dict=False)
     |      Returns tree in the form of a "bipartition dictionary"
     |
     |  branch_set(self)
     |      Returns set of (parent, child) tuples for all branches in tree
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
     |  build_remotechildren_dict(self)
     |      Constructs dict of all {parent:{remotechildren}} pairs in efficient manner.
     |      This dict can then be used directly or by remote_children() to speed up enquiries.
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
     |  cladedict(self, keep_remchild_dict=False)
     |      Returns tree in the form of a "clade dictionary"
     |
     |  cladegrep(self, pattern, minsize=2)
     |      Finds clades (monophyletic groups) where all leaves contain specified pattern
     |
     |  clear_attributes(self, remlist=None, keeplist=None)
     |
     |  clear_caches(self, preserve=[])
     |      Clears set of computed attributes. Use when tree structure is changed, and these
     |      are no longer reliable. Does not clear attributes listed in preserve
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
     |  copy_treeobject(self, copylengths=True, copyattr=True, interner=None)
     |      Returns copy of Tree object.
     |      copylengths: copy lengths (otherwise set to 0.0)
     |      copyattr: copy non-length attributes
     |      Caches are not copied.
     |      Similar to effect of copy.deepcopy but customized and much faster
     |
     |  deroot(self)
     |      If root is at bifurcation: remove root node, connect adjacent nodes
     |
     |  diameter(self, return_leaves=False)
     |      Return diameter: longest leaf-leaf distance along tree.
     |      If return_leaves is True: Return tuple with (maxdist, Leaf1, Leaf2)
     |
     |  find_bipart_nodes(self, leaves1, leaves2, check_arguments=True)
     |      Find two internal nodes corresponding to given bipartition of leaves.
     |
     |      Input: leaves1, leaves2: leaves on either side of internal branch (bipartition)
     |      Out: (basenode1, basenode2) nodes on either end of the branch corresponding to the bipartition
     |           basenode1 is at the base of leaves1 leaves
     |           basenode2 is at the base of leaves2 leaves
     |
     |      check_arguments: check that leaves1, leaves2 are valid bipartition of leaves
     |
     |      If tree is rooted at bifurcation, and bipartition corresponds to two halves of tree,
     |      then the two kids of root will be the basenodes (both basal branches are considered
     |      to be part of the same branch)
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
     |      Finds Most Recent Common Ancestor for the provided set of leaves.
     |      MRCA for a leaf is leaf itself
     |
     |  find_path_intersection(self, node1, node2)
     |      Returns intersection between paths from node1 and node2 to root (two-leaf mrca).
     |
     |  findbasenode(self, leafset)
     |      Finds node that is at the base of all leaves in leafset.
     |
     |  get_branch_attribute(self, node1, node2, attrname, default='')
     |      Get the value of any branch attribute.
     |      attrname: Name of attribute (e.g., "length")
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
     |  is_parent_child_pair(self, pnode, cnode)
     |      Returns True if cnode is child of pnode. False otherwise (also if order reversed)
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
     |  n_branches(self)
     |      Returns the number of branches in tree
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
     |  newick(self, printdist=True, printlabels=True, labelfield='label', precision=6, transdict=None, node_attributes=None, branch_attributes=None)
     |      Returns Newick format tree string representation of tree object, with optional metacomments
     |
     |  nexus(self, printdist=True, printlabels=True, labelfield='label', precision=6, translateblock=False, node_attributes=None, branch_attributes=None, colorlist=None, colorfg='#0000FF', colorbg='#000000')
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
     |  nodepath_to_branchset(self, nodepath)
     |      Input: nodepath (list of (node1, node2) tuples from .nodepath method)
     |      Output: set of branches (set of (parent, child) tuples)
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
     |  parsimony_assign_fits(self, fitpref=None)
     |      Performs parsimony analysis, and returns copy of tree object where one optimal fit
     |      (i.e., one giving the maximum parsimony number of changes) has been assigned to each node.
     |      fitpref = None: choose random fits at unlabelled nodes, where possible
     |      fitpref = <concrete state-string>: Set ambiguous state at unlabelled nodes to this value
     |               when possible (i.e., when compatible with MP score)
     |
     |  parsimony_count_changes(self, fitpref=None)
     |      Performs parsimony analysis on tree object, and returns parsimony score and
     |      dictionary giving number of changes in either direction:
     |          key = tuple (from,to), value = count in that direction
     |          {(state1, state2):count1_to_2,
     |           (state2, state1):count2_to_1}
     |
     |      fitpref = None: choose random fits at unlabelled nodes, where possible
     |      fitpref = <concrete state-string>: Set ambiguous state at unlabelled nodes to this value
     |               when possible (i.e., when compatible with MP score)
     |
     |  parsimony_possible_states(self)
     |      Performs parsimony analysis on tree object, and returns info about possible states
     |      at internal nodes.
     |
     |      Assumes tree object has .nodedict attribute: {nodeid:Nodestruct}, and that each
     |      Nodestruct has a .state attribute, which is either empty (unlabelled nodes)
     |      or a string specifying the state (e.g., "mink" or "human" or "A").
     |
     |      Returns: copy of tree object where nodestructs now also have attributes
     |      .primary_set, .secondary_set, and .optimal_set
     |      .optimal_set gives the full set of states that a node can have on the MP tree.
     |      Note: not all combinations of states will yield maximum parsimony score.
     |
     |      Based on Hartigan's 1973 parsimony algorithm.
     |      This version allows multifurcations, and internal nodes can have observed state.
     |
     |      Note: so far only implemented for one-character states (in the cladistic sense).
     |
     |  pathdist_as_ndarray(self, rooted=False)
     |      Return flattened version of pathdist_dict (in alphabetical leaf-pair order)
     |
     |  pathdist_dict(self, rooted=False)
     |      Returns dictionary giving path-distance (number of edges) for all pairs of leaves
     |
     |  patristic_distdict(self)
     |      Return nested dictionary giving all pairwise, patristic distances:
     |      dict[node1][node2] = patristic distance
     |
     |  possible_spr_prune_nodes(self)
     |      Utililty function when using spr function: where is it possible to prune
     |
     |  possible_spr_regraft_nodes(self, prune_node)
     |      Utility function when using spr function: where is it possible to regraft
     |      prune_node: the node below which pruning will take place (before regrafting)
     |
     |  prune_maxlen(self, nkeep, keeplist=[], return_leaves=False)
     |      Prune tree so remaining nkeep leaves spread out maximal percentage of branch length
     |      keeplist: optional list of leaves that must be included.
     |      Note: Best solution including keeplist may be less good than optimal solution
     |
     |  prune_subtree(self, basenode, keep_unary=False)
     |      Prune subtree rooted at basenode from self. Returns pruned subtree.
     |      keep_unary: keep left-behind, unary internal nodes in self (node with one descendant).
     |                  mostly useful when internal nodes have meaning beyond being hypothetical
     |                  ancestors (e.g., they can be observed, as in transmission trees)
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
     |      Removes named leaf from tree, cleans up so remaining tree structure is sane.
     |      Internal nodes may be removed, but will not change their IDs
     |      Note: assumes no unary nodes (leaf has at least one sibling)
     |
     |  remove_leaves(self, leaflist)
     |      Removes leaves in list from tree, cleans up so remaining tree structure is sane
     |
     |  rename_intnode(self, oldnum, newnum)
     |      Changes number of one internal node
     |
     |  rename_intnodes_to_match(self, other)
     |      Takes as input a tree (other) with the same topology as self, but with potentially
     |      different internal nodeIDs. Renames the internal nodeIDs in self so they are the
     |      same as those in other. Returns copy of self with new nodeIDs
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
     |  rootbip(self)
     |      For a tree rooted at a bifurcation: returns a tuple giving the following information
     |      about the bipartition on which the root is located:
     |              (Bipartition, leafset1, blen1, leafset2, blen2)
     |      where leafset1 and leafset2 are the two halves of the bipartition, and blen1 and blen2
     |      are the lengths of the branches leading from the root to their two basal nodes
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
     |  spr(self, prune_node=None, regraft_node=None)
     |      Subtree Pruning and Regrafting.
     |
     |      prune_node: basenode of subtree that will be pruned.
     |      regraft_node: node in remaining treestump below which subtree will be grafted
     |
     |      If no parameters are specified (both are None): perform random SPR
     |      If only prune_node is specified: choose random regraft_node
     |
     |      Must specify either both parameters, no parameters, or only prune_node
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
     |      Deprecated: Use treedist_RF instead.
     |      Compute symmetric tree distance (Robinson Foulds) between self and other tree.
     |      Normalised measure returned by default
     |
     |  treedist_RF(self, other, normalise=False, rooted=False, return_details=False)
     |      Compute symmetric tree distance (Robinson Foulds) between self and other tree.
     |      normalise: divide RF distance by the total number of bipartitions in the two trees
     |      rooted: take position of root into account.
     |      return_details: return list including intermediate values also:
     |          [treedist, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2,
     |           tree1_unique_biparts, tree2_unique_biparts, shared_biparts]
     |          where the last three are sets of bipartitions
     |
     |  treedist_pathdiff(self, other, rooted=False)
     |      Compute path difference tree-distance between self and other:
     |      Euclidean distance between nodepath-dist matrices considered as vectors.
     |      Measure described in M.A. Steel, D. Penny, Syst. Biol. 42 (1993) 126–141
     |      rooted=True: count traversal of root as 2 steps (not 1)
     |
     |  treesim(self, other, verbose=False)
     |      Compute normalised symmetric similarity between self and other tree
     |
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |
     |  from_biplist(biplist, interner=None)
     |      Constructor: Tree object from bipartition list
     |
     |  from_branchinfo(parentlist, childlist, lenlist=None, **attrlists)
     |      Constructor: Tree object from information about all branches in tree
     |
     |      Information about one branch is conceptually given as:
     |          parentnodeID, childnodeID, [length], [other named attributes...]
     |
     |      The function takes as input 2 or more separate lists containing:
     |          IDs of parents (internal nodes, so integer values)
     |          ID of children (internal or leaf nodes, so integer or string)
     |          Length of branches (optional)
     |          Any number of additional attributes given as <keyword=list>
     |
     |      The lists are assumed to have same length and be in same order (so index n in
     |      each list corresponds to same branch).
     |
     |      Note: most IDs appear multiple times in lists
     |      Note 2: can be used as workaround so user can specify IDs for internal nodes
     |
     |  from_cladedict(cladedict, interner=None)
     |
     |  from_leaves(leaflist, interner=None)
     |      Constructor: star-tree object from list of leaves
     |
     |  from_string(orig_treestring, transdict=None, interner=None, format='newick', label_attr_name='label', label_type=<class 'str'>)
     |      Constructor: Tree object from tree-string in specified format
     |
     |  from_topology(topology, interner=None)
     |      Constructor: Tree object from topology
     |
     |  randtree(leaflist=None, ntips=None, randomlen=False, name_prefix='s', interner=None)
     |      Constructor: tree with random topology from list of leaf names OR number of tips
     |
     |  ----------------------------------------------------------------------
     |  Static methods defined here:
     |
     |  is_property(cls, attr)
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  frozenset_leaves
     |
     |  leaf2index
     |
     |  nodedepthdict
     |
     |  parent_dict
     |      Lazy evaluation of _parent_dict when needed
     |
     |  remotechildren_dict
     |      Lazy evaluation of _remotechildren_dict when needed
     |
     |  rootdist
     |      Property (dictionary) giving the distance from each node in tree to the root
     |
     |  sorted_leaf_list
     |
     |  topology_bipart
     |      Returns set of Bipartitions representation of topology
     |
     |  topology_clade
     |      Returns set of Clades representation of topology
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

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
     |      list of weak references to the object
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
     |  __new__(*args, **kwargs) class method of builtins.Exception
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
     |  add_note(...)
     |      Exception.add_note(note) --
     |      add a note to the exception
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
     |  Class for storing and manipulating a number of trees, which all have the same leaves
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
     |  newick(self, printdist=True, printlabels=True, labelfield='label', precision=6, transdict=None, node_attributes=None, branch_attributes=None)
     |      Returns newick format tree for each tree in the set as a string
     |
     |  nexus(self, printdist=True, printlabels=True, labelfield='label', precision=6, translateblock=False, node_attributes=None, branch_attributes=None, colorlist=None, colorfg='#0000FF', colorbg='#000000')
     |      Returns nexus or format tree for each tree in the set
     |
     |  rootmid(self)
     |      Performs midpoint rooting on all trees in TreeSet
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  TreeSetIterator = <class 'phylotreelib.TreeSet.TreeSetIterator'>

    class TreeSummary(builtins.object)
     |  TreeSummary(trackbips=True, trackclades=False, trackroot=False, trackblen=False, trackdepth=False)
     |
     |  Class summarizing bipartitions and branch lengths (but not topologies) from many trees
     |
     |  Methods defined here:
     |
     |  __init__(self, trackbips=True, trackclades=False, trackroot=False, trackblen=False, trackdepth=False)
     |      TreeSummary constructor. Initializes relevant data structures
     |
     |  __len__(self)
     |
     |  add_tree(self, curtree, weight=1.0)
     |      Add tree object to treesummary, update all relevant bipartition summaries
     |
     |  compute_rootcred(self, tree)
     |      Returns root credibility (frequency of tree's root among observed trees) based
     |      on current root of sum_tree and information in self._rootbip_summary
     |
     |  contree(self, cutoff=0.5, allcompat=False)
     |      Returns a consensus tree built from selected bipartitions.
     |
     |  log_bipart_credibility(self, biptopology)
     |      Compute log bipartition credibility for topology (sum of log(freq) for all branches)
     |
     |  log_clade_credibility(self, cladetopology)
     |      Compute log clade credibility for topology (sum of log(freq) for all clades)
     |
     |  root_maxfreq(self, sum_tree)
     |      Uses info about root bipartitions in TreeSummary to place root on summary tree.
     |      Also sets tree attribute rootcred:
     |          probability (freq among input trees) of current location of root
     |
     |      If tree has branch lengths:
     |      Divides length of root bipartition among two branches in accordance with average
     |      fraction of lengths seen for this rootbip across all trees.
     |
     |  set_ca_node_depths(self, sum_tree, wt_count_burnin_filename_list)
     |      Set branch lengths on summary tree based on mean node depth for clades corresponding
     |      to MRCA of clade's leaves. (same as "--height ca" in BEAST's treeannotator)
     |      This means that all input trees are used when computing
     |      mean for each node (not just the input trees where that exact monophyletic clade
     |      is present)
     |
     |  set_clade_credibility(self, tree, precision=6)
     |      Set clade credibility on provided target tree based on freq of clade in TreeSummary.
     |
     |      NOTE: only works if all clades in tree have been observed at least once. The option
     |              will therefore not work with all rootings
     |
     |  set_mean_biplen(self, sum_tree)
     |      Sets branch-length, -var, and -sem for each branch in sum_tree,
     |      based on values in bipsummary.
     |      Should only be called when sum_tree constructed based on clades (MCC),
     |      but brlens are to be set based on mean bipartition values
     |
     |  set_mean_node_depths(self, sum_tree)
     |      Set branch lengths on summary tree based on mean node depth for clades corresponding
     |      to parent and child nodes (blen = depth_parent - depth_child).
     |
     |      NOTE 1: only meaningful if input trees are based on a clock model.
     |      NOTE 2: only works if all clades in tree have been observed at least once. The option
     |              will therefore not work with all rootings, and may also fail for majority rule
     |              consensus trees
     |      NOTE 3: only uses node depths from monophyletic clades (so some values may be set
     |      based on very few trees)
     |
     |  set_rootcredibility(self, sum_tree, precision=6)
     |      Returns sum_tree with root credibilities as attributes on each branch
     |      rootcred = fraction of trees in input set where the root was on this branch (bipartition)
     |      If root was never on a branch: assign the value 0.0
     |      Added as attribute .rootcred to Branchstruct for branches on sum_tree
     |      Also sets tree-attribute cumrootcred:
     |          sum of rootcredibilities for all branches (bipartitions) included on tree
     |
     |  update(self, other)
     |      Merge this object with external treesummary
     |
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |
     |  bipartsummary
     |      Property method for lazy evaluation of freq, var, sd, and sem for branches
     |
     |  cladesummary
     |      Property method for lazy evaluation of freq, var, and sem for node depths
     |
     |  rootbipsummary
     |      Property method for lazy evaluation of freq (=rootcred) for rootbips
     |
     |  sorted_biplist
     |      Return list of bipartitions.
     |      First external (leaf) bipartitions sorted by leafname.
     |      Then internal bipartitions sorted by freq
     |
     |  sorted_rootbips
     |      Return list of root-bipartitions (branches where root has been seen), sorted by
     |      occurrence (count on tree samples added to TreeSummary)
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

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
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

    class TreefileBase(builtins.object)
     |  TreefileBase(filename, filecontent, interner)
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
     |  __init__(self, filename, filecontent, interner)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  close(self)
     |      For explicit closing of Treefile before content exhausted
     |
     |  get_treestring(self)
     |      Return next tree-string
     |
     |  readtree(self, returntree=True)
     |      Reads one treestring from file and returns as Tree object if requested.
     |      Returns None when exhausted file
     |
     |  readtrees(self, discardprop=0.0)
     |      Reads trees from file and returns as TreeSet object. Can discard fraction of trees
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables
     |
     |  __weakref__
     |      list of weak references to the object

FUNCTIONS
    count_bytestring(filename, bytestring)
        Fast counting of specific pattern. Bytestring argument must be given
        with b modifier (e.g., b');')

    count_trees_by_parsing(filename, args)

    fast_treecount(filename, fileformat='nexus')
        Heuristic: count patterns ([;=
        ] etc) to infer number of trees

    main()
        # Placeholder: Insert test code here and run module in standalone mode

    remove_comments(text)
        Takes input string and strips away commented text, delimited by '[' and ']'.
        Also deals with nested comments.

FILE
    /Users/agpe/Documents/code/python/phylotreelib/phylotreelib.py



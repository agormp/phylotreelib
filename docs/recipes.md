## phylotreelib recipes

This file contains task-oriented examples for common `phylotreelib` workflows.

The recipes are intentionally practical and somewhat redundant. The goal is to provide blocks of code that can be copied and easily adapted for related tasks.

### Terminology note: tree orientation and "height" vs "depth"

**Note for users upgrading from version 1:** In **phylotreelib** v. 2, trees are now treated as *having the root at the top and leaves at the bottom* — matching the standard convention in most phylogeny software and computer science.

Accordingly, it uses:

- **height(node)** = distance from the **tips (leaves)** to the node (i.e., "time before most recent leaf")

In v1.x, this concept was called **depth**. All depth-based names were renamed to height-based names in v2.0.0. See [Upgrading from 1.x to 2.x](#upgrading-from-1x-to-2x) for the full rename table.

---

### 1. Reading trees

#### Read a single tree from a NEXUS file

```python
import phylotreelib as pt
with pt.Nexustreefile("example.nexus") as tf:
    tree = tf.readtree()
```

#### Read a single tree from a Newick file

```python
import phylotreelib as pt
with pt.Newicktreefile("example.newick") as tf:
    tree = tf.readtree()
```

#### Let phylotreelib autodetect the tree file format

```python
import phylotreelib as pt
with pt.Treefile("example.tree") as tf:
    tree = tf.readtree()
```

#### Iterate over all trees in a file without storing them all in memory

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/primate-mtDNA.trees"
urllib.request.urlretrieve(URL, "primate-mtDNA.trees")

import phylotreelib as pt
with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        # do something with one tree at a time
        print(f"n_tips={len(tree.leaves):>3d}  length={tree.length():.6g}")
        break  # remove if you want to iterate over all trees
```

#### Read all trees into a TreeSet

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/primate-mtDNA.trees"
urllib.request.urlretrieve(URL, "primate-mtDNA.trees")

import phylotreelib as pt
with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    treeset = tf.readtrees()
print(len(treeset))
```

---

### 2. Constructing trees directly

#### From a Newick string

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
```

#### From a list of leaf names (star tree)

```python
import phylotreelib as pt
tree = pt.Tree.from_leaves(["A", "B", "C", "D"])
```

#### From explicit branch information

```python
import phylotreelib as pt
parents = [100, 100, 101, 101, 102, 102]
children = [101, 102, "A", "B", "C", "D"]
lengths = [1.0, 1.0, 2.0, 3.0, 2.0, 3.0]
labels = [0.95, 0.65, 0.99, 0.78, 0.93, 0.99]
tree = pt.Tree.from_branchinfo(
    parentlist=parents,
    childlist=children,
    lenlist=lengths,
    label=labels,
)
print(tree)
```

Output:

```text
+--------------------------------------------+
|  parent  |  child  |  branchlen  |  label  |
+--------------------------------------------+
|     100  |    101  |          1  |   0.95  |
|     100  |    102  |          1  |   0.65  |
|     101  |      A  |          2  |   0.99  |
|     101  |      B  |          3  |   0.78  |
|     102  |      C  |          2  |   0.93  |
|     102  |      D  |          3  |   0.99  |
+--------------------------------------------+

4 Leaves:
-
A
B
C
D
```

#### Random tree

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=20, randomlen=True)
```

#### Random tree with specified leaf names

```python
import phylotreelib as pt
tree = pt.Tree.randtree(leaflist=["Human", "Chimp", "Bonobo", "Gorilla", "Orangutan"], randomlen=True)
```


---

### 3. Querying tree structure

#### Find the MRCA of a set of leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
mrca = tree.find_mrca({"A", "B"})
print(mrca)
```

#### Get all descendant leaves ("remote children") below a node

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,((C:1,D:1):1,E:2):1);")
mrca = tree.find_mrca({"C", "E"})
print(tree.remote_children(mrca))
```

Note: .remote_children() returns a set of leaf names.

#### Get immediate descendant nodes ("children") below a node

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,((C:1,D:1):1,E:2):1);")
mrca = tree.find_mrca({"C", "E"})
print(tree)
print(f"Internal node: {mrca}")
print(f"Children of internal node: {list(tree.children(mrca))}")
```

Note: .children() returns a dict_keys object which behaves as a set (you can iterate over it and perform set arithmetic)

#### Find leaves near a focal leaf by patristic distance

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
neighbors = tree.nearleafs("A", maxdist=2.1)
print(neighbors)
```

#### Find the nearest N leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
print(tree.nearest_n_leaves("A", 3))
```

#### Compute patristic distance between two nodes

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
print(tree.nodedist("A", "D"))
```

#### Get the path along the tree between two nodes

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
path = tree.nodepath("A", "D")
print(path)
```

#### Find diameter and two farthest leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
maxdist, leaf1, leaf2 = tree.diameter(return_leaves=True)
print(maxdist, leaf1, leaf2)
```

#### Find a central or common leaf

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
leafset = ["A", "B", "C", "D"]
print(tree.find_central_leaf(leafset))
print(tree.find_common_leaf(leafset))
```

---

### 4. Rooting and rerooting

#### Midpoint rooting

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5, randomlen=True)
print(tree)
tree.rootmid()
print(tree)
```

#### Minimum-variance rooting

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5, randomlen=True)
tree.rootminvar()
```

#### Rooting on an outgroup with one taxon

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(Gorilla:1.5,((Human:2,(Chimpanzee:1,Bonobo:1):1):1,Orangutan:5):1.5);")
print(tree)
tree.rootout("Orangutan")
print(tree)
```

Note: Outgroup rooting places the root on the branch leading to the outgroup. The original basal branch length is split evenly, so half is assigned to the outgroup branch and half to the ingroup branch.

#### Rooting on an outgroup with multiple taxa

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(Gorilla:1.5,((Human:2,(Chimpanzee:1,Bonobo:1):1):1,Orangutan:5):1.5);")
og = set(["Chimpanzee","Bonobo"])
tree.rootout(og)
print(tree)
```

#### Reroot on a specified branch

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
node = tree.find_mrca({"C", "D"})
parent = tree.parent(node)

print(tree)
print(f"Add new root node on edge {parent}->{node} at 25% from {parent}")
tree.reroot(node1=parent, node2=node, frac_from_node1=0.25)
print(tree)
```

#### Remove one leaf

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
print(tree)
tree.remove_leaf("D")
print(tree.newick())
```

#### Remove several leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);")
tree.remove_leaves(["B", "F"])
print(tree.newick())
```

#### Extract a subtree

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
print(tree)

# `subtree()` returns (subtree, basalbranch)
# - subtree is a Tree when basenode is internal
# - subtree is a leaf name (str) when basenode is a leaf
sub, basal = tree.subtree(tree.find_mrca({"C", "D"}))
print("Basal branch length:", basal.length)

# Print the extracted subtree (internal-node case)
print(sub.newick())
```

**Note:** `subtree()` does not modify the original tree. It returns a detached copy of the subtree
and a copy of the incoming `Branchstruct` (the basal branch into the basenode).

#### Prune a subtree out of the tree

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
print(tree)

# `prune_subtree()` returns (subtree, basalbranch) and removes the subtree from `tree`
sub, basal = tree.prune_subtree(tree.find_mrca({"C", "D"}))
print("Basal branch length that was removed:", basal.length)

print("Remaining tree:")
print(tree)

print("Pruned subtree:")
if isinstance(sub, str):
    print(sub)          # pruned a single leaf
else:
    print(sub.newick()) # pruned an internal subtree
```

#### Insert an internal node (split off descendants)

This operation *splits off* one or more existing child branches under a fresh internal node.
Only the new upper branch gets the provided Branchstruct; the moved child branches keep their
original Branchstructs (including their lengths).

```python
import phylotreelib as pt
from phylotreelib import Branchstruct

tree = pt.Tree.from_string("(A:1,B:1,C:1,D:1);")
print(tree)

# Move A and B under a new internal node below the root (node 0)
newnode = tree.split_off_children(0, ["A", "B"], Branchstruct(length=0.5, label="newsplit"))
print("New internal node:", newnode)
print(tree)
```

**Note:** `insert_node()` is a deprecated alias for `split_off_children()`.

#### Insert a node on an existing branch

This operation inserts a node *on* a specific existing edge (parent → child) and splits the
original branch length between the two new segments.

```python
import phylotreelib as pt

tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
print(tree)

node = tree.find_mrca({"C", "D"})
parent = tree.parent(node)

# Insert a node halfway along the incoming branch to `node`
newnode = tree.add_node_on_branch(parent, node, frac_from_parent=0.5)
print("Inserted node:", newnode)
print(tree)
```



#### Graft one tree onto another

`graft()` inserts a graft point on an existing edge in the target tree (parent → child), then
attaches either a whole Tree (via its root) or a single leaf name below that graft point.

```python
import phylotreelib as pt

left = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
right = pt.Tree.from_string("(X:1,(Y:1,Z:1):1);")

# Choose the edge in `left` where the graft point should be inserted:
child = left.find_mrca({"C", "D"})
parent = left.parent(child)

print(left)

# Attach `right` below the graft point, 30% of the way from parent to child
# The new connecting branch above right subtree has length 0.25
left.graft(right, parent=parent, child=child,
           frac_from_parent=0.3, connect_length=0.25,
           graftlabel="grafted_")

print(f"right grafted between nodes {parent} and {child}")
print(left)
```

#### SPR with specified prune and regraft points

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
prune_node = tree.find_mrca({"C", "D"})
regraft_node = tree.find_mrca({"E", "F"})
print(tree.newick(printdist=False))
tree.spr(prune_node=prune_node, regraft_node=regraft_node)
print(tree.newick(printdist=False))
```

#### SPR with specified prune point and random regraft point

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5, randomlen=True)
print(tree.newick(precision=2))
prune_node = next(iter(tree.leaves))
prune_node, regraft_node = tree.spr(prune_node=prune_node)
print(f"prune_node; {prune_node}\tregraft_node = {regraft_node}")
print(tree.newick(precision=2))
```

#### Fully random SPR moves

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=7, randomlen=True)
print(f"{tree.newick(printdist=False)}\tlength: {tree.length()}")
for _ in range(10):
    tree.spr()
    print(f"{tree.newick(printdist=False)}\tlength: {tree.length()}")
```

---

### 6. Representative pruning

#### Keep leaves that maximize retained tree length

```python
# Download a large example tree (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/large_njtree.txt"
urllib.request.urlretrieve(URL, "large_njtree.txt")

import phylotreelib as pt

# large_njtree.txt contains a single Newick tree (many tips)
with open("large_njtree.txt") as fh:
    tree = pt.Tree.from_string(fh.read().strip())

tree.prune_maxlen(nkeep=50)
print(f"n_leaves={len(tree.leaves)}  length={tree.length():.6g}")
```

#### Keep leaves that maximize retained tree length while forcing some leaves to remain

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=200, randomlen=True)
tree.prune_maxlen(nkeep=30, keeplist=["s001", "s120", "s178"])
```

#### Ask which leaves would be kept without actually pruning, color those on original tree

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=50, randomlen=True)
tree.rootminvar()
keepers = tree.prune_maxlen(nkeep=10, return_leaves=True)
with open("colorprune.nexus", "w") as outfile:
    outfile.write(tree.nexus(colorlist=keepers, colorfg="#FF0000"))
```

Note: the leaves selected to be retained after pruning will now be colored red when viewed in FigTree (and possibly other viewers)

---

### 7. Comparing trees

#### Robinson–Foulds distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=30, randomlen=True)
t2 = t1.copy_treeobject()
t2.spr()
print(t1.treedist_RF(t2))
print(t1.treedist_RF(t2, normalise=True))
```

#### Rooted RF distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=20, randomlen=True)
t2 = t1.copy_treeobject()
t2.rootmid()
print(t1.treedist_RF(t2, rooted=False))
print(t1.treedist_RF(t2, rooted=True))
```

#### Path-difference distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=30, randomlen=True)
t2 = t1.copy_treeobject()
for _ in range(3):
    t2.spr()
print(t1.treedist_pathdiff(t2))
```

#### Request detailed RF output

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=10, randomlen=True)
t2 = t1.copy_treeobject()
t2.spr()
(treedist, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2, tree1_unique_biparts,
 tree2_unique_biparts, shared_biparts) = t1.treedist_RF(t2, return_details=True)
print(treedist, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2)
```

---

### 8. Distance matrices

#### Build neighbor-joining and UPGMA trees from a nested distance dictionary

```python
import phylotreelib as pt
distdict = {
    "A": {"A": 0, "B": 5, "C": 9, "D": 9},
    "B": {"A": 5, "B": 0, "C": 10, "D": 10},
    "C": {"A": 9, "B": 10, "C": 0, "D": 8},
    "D": {"A": 9, "B": 10, "C": 8, "D": 0},
}
dmat = pt.Distmatrix.from_distdict(distdict)
print(dmat.nj().newick())
print(dmat.upgma().newick())
```

#### Build neighbor-joining tree from a NumPy array

```python
import numpy as np
import phylotreelib as pt
arr = np.array([
    [0, 5, 9, 9],
    [5, 0, 10, 10],
    [9, 10, 0, 8],
    [9, 10, 8, 0],
], dtype=float)
names = ["A", "B", "C", "D"]
dmat = pt.Distmatrix.from_numpy_array(arr, names)
print(dmat.nj().newick())
```

---

### 9. Consensus trees and other tree summaries

**Note:** For command-line use, the separate tool `sumt` provides the same functionality with a simpler interface. The examples below show the explicit low-level workflow in `phylotreelib` (TreeSummary --> SummaryTreeBuilder --> TreePostProcessor --> PrintSpec).

* `sumt` GitHub page: [https://github.com/agormp/sumt](https://github.com/agormp/sumt/)

More detailed discussion of summary tree methods is in `docs/api.md`.

#### Build a majority-rule consensus tree

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/primate-mtDNA.trees"
urllib.request.urlretrieve(URL, "primate-mtDNA.trees")

import phylotreelib as pt

# 1) Accumulate bipartition frequencies + mean branch lengths
ts = pt.TreeSummary(trackbips=True, trackblen=True)

with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

# 2) Build majority-rule consensus topology
stb = pt.SummaryTreeBuilder(ts)
contree = stb.contree(cutoff=0.5, allcompat=False)

# 3) Annotate (support + branch-length stats)
tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(contree)

# 4) Printing (optional)
pt.configure_sumtree_printing(contree, treetype="con", blen="biplen", print_meta=True, precision=7)

print(contree.newick())
```

#### Build a majority-rule consensus tree with all compatible bipartitions added

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/primate-mtDNA.trees"
urllib.request.urlretrieve(URL, "primate-mtDNA.trees")

import phylotreelib as pt

ts = pt.TreeSummary(trackbips=True, trackblen=True)

with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

stb = pt.SummaryTreeBuilder(ts)
contree = stb.contree(cutoff=0.5, allcompat=True)

tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(contree)

pt.configure_sumtree_printing(contree, treetype="all", blen="biplen", print_meta=True, precision=7)

print(contree.newick())
```

#### Build an MCC tree with common-ancestor heights

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/primate-mtDNA.trees"
urllib.request.urlretrieve(URL, "primate-mtDNA.trees")

import phylotreelib as pt

# Pass 1: collect clade frequencies + observed topologies (needed for MCC selection)
ts = pt.TreeSummary(trackclades=True, tracktopo=True)

with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

stb = pt.SummaryTreeBuilder(ts)
mcctree = stb.max_clade_cred_tree()

# Pass 2: compute common-ancestor heights (CA-height) for nodes on the chosen topology
plan = pt.CAHeightEstimator.build_plan(mcctree, trackci=False)
est = pt.CAHeightEstimator(plan, trackci=False)

with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        est.add_tree(tree)

est.write_into(mcctree)
mcctree.set_blens_from_heights()

# Annotate support etc from TreeSummary (does not change heights you just wrote)
tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(mcctree)

pt.configure_sumtree_printing(mcctree, treetype="mcc", blen="caheight", print_meta=True, precision=7)

print(mcctree.newick())
```

#### Build an MBC tree and mid-point root it

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/primate-mtDNA.trees"
urllib.request.urlretrieve(URL, "primate-mtDNA.trees")

import phylotreelib as pt

# Need to keep track of bipartitions + observed topologies to construct MBC
ts = pt.TreeSummary(trackbips=True, tracktopo=True)

with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

stb = pt.SummaryTreeBuilder(ts)
mbctree = stb.max_bipart_cred_tree()

# MBC is unrooted; choose a root for presentation if desired
mbctree.rootmid()

tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(mbctree)

pt.configure_sumtree_printing(mbctree, treetype="mbc", blen="input", print_meta=True, precision=7)

print(mbctree.newick())
```


---

### 10. Exporting with labels, metadata, and color

#### Plain Newick output

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
with open("tree.newick", "w") as out:
    out.write(tree.newick())
```

#### NEXUS output with translate block

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
with open("tree.nexus", "w") as out:
    out.write(tree.nexus(translateblock=True))
```

File contents:

```
#NEXUS

begin trees;
    translate
        1     A,
        2     B,
        3     C,
        4     D
    ;
	tree nexus_tree = ((1:1,2:1):1,(3:1,4:1):1);
end;
```

#### Print branch metadata in BEAST/FigTree-style comments

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
node = tree.find_mrca({"A", "B"})
parent = tree.parent(node)
tree.set_branch_attribute(parent, node, "posterior", 0.98)
pt.configure_basic_printing(
    tree,
    print_meta=True,
    branch_attrs=("posterior",),
    precision=6,
)
print(tree.newick())
```

Output:

```
((A:1,B:1):1[&posterior=0.98],(C:1,D:1):1);
```

#### Color leaves selected by a clade query

```python
import phylotreelib as pt
tree = pt.Tree.from_string(
    "((Mink_1:1,Mink_2:1):1,(Human_1:1,Human_2:1):1,(Bat_1:1,Bat_2:1):1);"
)
mrca = tree.find_mrca({"Mink_1", "Mink_2"})
red_leaves = sorted(tree.remote_children(mrca))
with open("colored.nexus", "w") as out:
    out.write(
        tree.nexus(
            colorlist=red_leaves,
            colorfg="#FF0000",
            colorbg="#000000",
        )
    )
```

#### Color leaves selected by name pattern

```python
import phylotreelib as pt
tree = pt.Tree.from_string(
    "((Mink_1:1,Mink_2:1):1,(Human_1:1,Human_2:1):1,(Bat_1:1,Bat_2:1):1);"
)
red_leaves = [leaf for leaf in tree.leaves if leaf.startswith("Human_")]
with open("colored_humans.nexus", "w") as out:
    out.write(tree.nexus(colorlist=red_leaves, colorfg="#0000FF"))
```

---

### 11. Error handling

```python
import phylotreelib as pt
try:
    with pt.Newicktreefile("broken_treefile.newick") as tf:
        tree = tf.readtree()
except pt.TreeError as err:
    print("Tree-related error:", err)
```

---

# phylotreelib recipes

This file contains task-oriented examples for common `phylotreelib` workflows.

The recipes are intentionally practical and somewhat redundant. The goal is that you can copy a block, adapt filenames or leaf names, and get work done quickly.

---

# 1. Reading trees

## Read a single tree from a NEXUS file

```python
import phylotreelib as pt
with pt.Nexustreefile("example.nexus") as tf:
    tree = tf.readtree()
```

## Read a single tree from a Newick file

```python
import phylotreelib as pt
with pt.Newicktreefile("example.newick") as tf:
    tree = tf.readtree()
```

## Let phylotreelib autodetect the format

```python
import phylotreelib as pt
with pt.Treefile("example.tree") as tf:
    tree = tf.readtree()
```

## Iterate over all trees in a file without storing them all in memory

```python
import phylotreelib as pt
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        # do something with one tree at a time
        pass
```

## Read all trees into a TreeSet

```python
import phylotreelib as pt
with pt.Nexustreefile("posterior.trees") as tf:
    treeset = tf.readtrees()
print(len(treeset))
```

---

# 2. Constructing trees directly

## From a Newick string

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
```

## From a list of leaf names (star tree)

```python
import phylotreelib as pt
tree = pt.Tree.from_leaves(["A", "B", "C", "D"])
```

## From explicit branch information

```python
import phylotreelib as pt
parents = [100, 100, 101, 101, 102, 102]
children = [101, 102, "A", "B", "C", "D"]
lengths = [1.0, 1.0, 2.0, 3.0, 2.0, 3.0]
labels = ["", "", "", "", "", ""]
tree = pt.Tree.from_branchinfo(
    parentlist=parents,
    childlist=children,
    lenlist=lengths,
    label=labels,
)
print(tree)
```

## Random tree for testing

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=20, randomlen=True)
```

---

# 3. Querying tree structure

## Find the MRCA of a set of leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
mrca = tree.find_mrca({"A", "B"})
print(mrca)
```

## Get all descendant leaves below a node

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
mrca = tree.find_mrca({"C", "D"})
print(tree.remote_children(mrca))
```

## Find leaves near a focal leaf by patristic distance

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
neighbors = tree.nearleafs("A", maxdist=2.1)
print(neighbors)
```

## Find the nearest N leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
print(tree.nearest_n_leaves("A", 3))
```

## Compute pairwise distance between two leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
print(tree.nodedist("A", "D"))
```

## Get the path between two nodes

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
path = tree.nodepath("A", "D")
print(path)
```

## Find diameter and the two farthest leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
maxdist, leaf1, leaf2 = tree.diameter(return_leaves=True)
print(maxdist, leaf1, leaf2)
```

## Find a central or common leaf

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
leafset = ["A", "B", "C", "D"]
print(tree.find_central_leaf(leafset))
print(tree.find_common_leaf(leafset))
```

---

# 4. Rooting and rerooting

## Midpoint rooting

```python
import phylotreelib as pt
with pt.Newicktreefile("unrooted.newick") as tf:
    tree = tf.readtree()
tree.rootmid()
```

## Minimum-variance rooting

```python
import phylotreelib as pt
with pt.Newicktreefile("unrooted.newick") as tf:
    tree = tf.readtree()
tree.rootminvar()
```

## Rooting on an outgroup

```python
import phylotreelib as pt
with pt.Newicktreefile("unrooted.newick") as tf:
    tree = tf.readtree()
tree.rootout({"Out1", "Out2"})
```

## Reroot on a specified branch

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
node = tree.find_mrca({"C", "D"})
parent = tree.parent(node)
tree.reroot(parent, node, node1dist=0.4)
```

---

# 5. Modifying topology

## Remove one leaf

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
tree.remove_leaf("D")
```

## Remove several leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);")
tree.remove_leaves(["E", "F"])
```

## Extract a subtree

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
basenode = tree.find_mrca({"A", "B", "C", "D"})
sub = tree.subtree(basenode)
print(sub.newick())
```

## Prune a subtree out of the tree

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
basenode = tree.find_mrca({"E", "F"})
pruned = tree.prune_subtree(basenode)
print("Remaining tree:", tree.newick())
print("Pruned subtree:", pruned.newick())
```

## Insert an internal node

```python
import phylotreelib as pt
from phylotreelib import Branchstruct
tree = pt.Tree.from_string("(A:1,B:1,C:1,D:1);")
newnode = tree.insert_node(0, ["A", "B"], Branchstruct(length=0.5, label="newsplit"))
print(newnode)
print(tree.newick())
```

## Graft one tree onto another

```python
import phylotreelib as pt
left = pt.Tree.from_string("((A:1,B:1):1,C:1);")
right = pt.Tree.from_string("(X:1,Y:1);")
node = left.find_mrca({"A", "B"})
left.graft(right, node1=node, blen1=0.2, blen2=0.3, graftlabel="g2")
print(left.newick())
```

## SPR with specified prune and regraft points

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
prune_node = tree.find_mrca({"C", "D"})
regraft_node = tree.find_mrca({"E", "F"})
tree.spr(prune_node=prune_node, regraft_node=regraft_node)
print(tree.newick())
```

## SPR with specified prune point and random regraft point

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=10, randomlen=True)
prune_node = next(iter(tree.intnodes))
tree.spr(prune_node=prune_node)
```

## Fully random SPR moves

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=20, randomlen=True)
for _ in range(10):
    tree.spr()
```

---

# 6. Representative pruning

## Keep leaves that maximize retained tree length

```python
import phylotreelib as pt
with pt.Nexustreefile("big.tree") as tf:
    tree = tf.readtree()
small = tree.prune_maxlen(nkeep=50)
```

## Keep leaves that maximize retained tree length while forcing some leaves to remain

```python
import phylotreelib as pt
with pt.Nexustreefile("big.tree") as tf:
    tree = tf.readtree()
small = tree.prune_maxlen(
    nkeep=50,
    keeplist=["human1", "human2", "mink1", "mink2", "mink3"],
)
```

## Keep approximately evenly spaced leaves

```python
import phylotreelib as pt
with pt.Nexustreefile("big.tree") as tf:
    tree = tf.readtree()
small = tree.numberprune(nkeep=50)
```

## Ask which leaves would be kept without actually pruning

```python
import phylotreelib as pt
with pt.Nexustreefile("big.tree") as tf:
    tree = tf.readtree()
keepers = tree.prune_maxlen(nkeep=50, return_leaves=True)
print(sorted(keepers))
```

---

# 7. Comparing trees

## Robinson–Foulds distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=30, randomlen=True)
t2 = t1.copy_treeobject()
for _ in range(5):
    t2.spr()
print(t1.treedist_RF(t2))
print(t1.treedist_RF(t2, normalise=True))
```

## Rooted RF distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=20, randomlen=True)
t2 = t1.copy_treeobject()
t2.rootmid()
print(t1.treedist_RF(t2, rooted=True))
```

## Path-difference distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=30, randomlen=True)
t2 = t1.copy_treeobject()
for _ in range(3):
    t2.spr()
print(t1.treedist_pathdiff(t2))
```

## Request detailed RF output

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=10, randomlen=True)
t2 = t1.copy_treeobject()
t2.spr()
details = t1.treedist_RF(t2, return_details=True)
print(details)
```

---

# 8. Distance matrices

## Build a tree from a nested distance dictionary

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

## Build a tree from a NumPy array

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

# 9. Posterior tree summaries

## Build a majority-rule consensus tree

```python
import phylotreelib as pt
ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)
contree = ts.contree(cutoff=0.5)
print(contree.newick())
```

## Build a majority-rule consensus tree with all compatible bipartitions added

```python
import phylotreelib as pt
ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)
contree = ts.contree(cutoff=0.5, allcompat=True)
print(contree.newick())
```

## Build an MCC tree

```python
import phylotreelib as pt
ts = pt.TreeSummary(trackclades=True, tracktopo=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)
mcctree = ts.max_clade_cred_tree()
print(mcctree.newick())
```

## Build an MBC tree

```python
import phylotreelib as pt
ts = pt.TreeSummary(trackbips=True, tracktopo=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)
mbctree = ts.max_bipart_cred_tree()
print(mbctree.newick())
```

## Build a summary tree with one convenience call

```python
import phylotreelib as pt
ts = pt.TreeSummary(
    trackbips=True,
    trackclades=True,
    trackroot=True,
    trackblen=True,
    trackdepth=True,
    trackrootblen=True,
    tracktopo=True,
)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)
sumtree = ts.compute_sumtree(
    treetype="mcc",
    blen="depth",
    rooting="maxfreq",
)
```

## Use the top-level helper instead of the bound method

```python
import phylotreelib as pt
ts = pt.TreeSummary(
    trackbips=True,
    trackclades=True,
    trackroot=True,
    trackblen=True,
    trackdepth=True,
    trackrootblen=True,
    tracktopo=True,
)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)
sumtree = pt.build_sumtree(ts, treetype="con", blen="biplen", rooting="mid")
```

---

# 10. Exporting with labels, metadata, and color

## Plain Newick output

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
with open("tree.newick", "w") as out:
    out.write(tree.newick())
```

## NEXUS output with translate block

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
with open("tree.nexus", "w") as out:
    out.write(tree.nexus(translateblock=True))
```

## Print branch metadata in BEAST/FigTree-style comments

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

## Change which branch attribute is printed as the internal label

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
node = tree.find_mrca({"A", "B"})
parent = tree.parent(node)
tree.set_branch_attribute(parent, node, "bipartition_cred", 0.95)
tree.set_print_spec(labelfield="bipartition_cred", printlabels=True)
print(tree.newick())
```

## Color leaves selected by a clade query

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

## Color leaves selected by name pattern

```python
import phylotreelib as pt
tree = pt.Tree.from_string(
    "((Mink_1:1,Mink_2:1):1,(Human_1:1,Human_2:1):1,(Bat_1:1,Bat_2:1):1);"
)
red_leaves = sorted([leaf for leaf in tree.leaves if leaf.startswith("Human_")])
with open("colored_humans.nexus", "w") as out:
    out.write(tree.nexus(colorlist=red_leaves, colorfg="#0000FF"))
```

---

# 11. Error handling

```python
import phylotreelib as pt
try:
    with pt.Newicktreefile("broken_treefile.newick") as tf:
        tree = tf.readtree()
except pt.TreeError as err:
    print("Tree-related error:", err)
```

---

# 12. When to use which summary-tree method

This is just a quick practical guide.

- Use **consensus trees** when you want a summary of commonly observed splits.
- Use **all-compatible consensus** when you want a more resolved consensus tree while still respecting compatibility.
- Use **MCC trees** when you want the sampled tree whose clades have the highest joint support.
- Use **MBC trees** when branch/bipartition support is the primary criterion.
- Use **HIPSTR** when you want a summary tree based on independent posterior subtree reconstruction.

More detailed discussion is in `docs/api.md`.


## 9. Consensus trees and other tree summaries

If you mainly want “TreeAnnotator-like” functionality from posterior tree samples, the command-line
tool **`sumt`** is usually the fastest and simplest way to get a consensus/MCC/MBC/HIPSTR tree,
with branch-length summaries and annotations already configured.

- `sumt` (CLI) on GitHub: https://github.com/agormp/sumt

The examples below show the *explicit* (low-level) workflow inside `phylotreelib`. Nothing is
hidden: you can see exactly which step constructs the topology and which step adds annotations.

### 9.1 Build a majority-rule consensus tree (bipartition support)

```python
import phylotreelib as pt

# 1) Accumulate bipartitions (+ optional branch-length stats)
ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

# 2) Topology from bipartitions (majority rule)
stb = pt.SummaryTreeBuilder(ts)
contree = stb.contree(cutoff=0.5, allcompat=False)

# 3) Paint/annotate the tree with tracked attributes
tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(contree)   # sets bipartition_cred (+ length/length_sd if tracked)

# 4) Printing: support as internal labels, optional meta-comments
contree.set_print_spec(
    labelfield="bipartition_cred",
    printlabels=True,
    printdist=True,
    print_meta=False,
    precision=7,
)

print(contree.newick())
```

### 9.2 Consensus tree with all compatible bipartitions added

```python
import phylotreelib as pt

ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

stb = pt.SummaryTreeBuilder(ts)
contree = stb.contree(cutoff=0.5, allcompat=True)

tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(contree)

contree.set_print_spec(labelfield="bipartition_cred", precision=7)
print(contree.newick())
```

### 9.3 MCC tree with “common ancestor depths” (TreeAnnotator --height ca)

This requires a **second pass over the tree samples** to compute common-ancestor depths (CA depths).

```python
import phylotreelib as pt

# 1) Track clades + topologies (for MCC)
ts = pt.TreeSummary(trackclades=True, tracktopo=True)

with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

# 2) Build MCC topology
stb = pt.SummaryTreeBuilder(ts)
mcctree = stb.max_clade_cred_tree()

# 3) CA depths via CADepthEstimator (stream trees again)
plan = pt.CADepthEstimator.build_plan(mcctree, trackci=False)
est = pt.CADepthEstimator(plan, trackci=False)

with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        est.add_tree(tree)

est.write_into(mcctree)         # writes depth/depth_sd onto Nodestructs
mcctree.set_blens_from_depths()

# 4) Annotate supports etc.
tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(mcctree)

mcctree.set_print_spec(
    print_meta=True,
    node_attrs=("clade_cred", "depth", "depth_sd"),
    branch_attrs=("length",),
    precision=7,
)
print(mcctree.newick())
```

### 9.4 MBC tree and midpoint rooting

```python
import phylotreelib as pt

# Need bipartitions + topologies to choose the maximum-bipartition-credibility topology
ts = pt.TreeSummary(trackbips=True, tracktopo=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

stb = pt.SummaryTreeBuilder(ts)
mbctree = stb.max_bipart_cred_tree()

# Root the chosen topology (Tree methods)
mbctree.rootmid()

tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(mbctree)

mbctree.set_print_spec(labelfield="bipartition_cred", precision=7)
print(mbctree.newick())
```

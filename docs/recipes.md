## phylotreelib recipes

This file contains task-oriented examples for common `phylotreelib` workflows.

The recipes are intentionally practical and somewhat redundant. The goal is to provide blocks of code that can be copied and easily adapted for related tasks.

---

## 1. Reading trees

### Read a single tree from a NEXUS file

```python
import phylotreelib as pt
with pt.Nexustreefile("example.nexus") as tf:
    tree = tf.readtree()
```

### Read a single tree from a Newick file

```python
import phylotreelib as pt
with pt.Newicktreefile("example.newick") as tf:
    tree = tf.readtree()
```

### Let phylotreelib autodetect the tree file format

```python
import phylotreelib as pt
with pt.Treefile("example.tree") as tf:
    tree = tf.readtree()
```

### Iterate over all trees in a file without storing them all in memory

```python
import phylotreelib as pt
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        # do something with one tree at a time
        pass
```

### Read all trees into a TreeSet

```python
import phylotreelib as pt
with pt.Nexustreefile("posterior.trees") as tf:
    treeset = tf.readtrees()
print(len(treeset))
```

---

## 2. Constructing trees directly

### From a Newick string

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
```

### From a list of leaf names (star tree)

```python
import phylotreelib as pt
tree = pt.Tree.from_leaves(["A", "B", "C", "D"])
```

### From explicit branch information

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

### Random tree

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=20, randomlen=True)
```

### Random tree with specified leaf names

```python
import phylotreelib as pt
tree = pt.Tree.randtree(leaflist=["Human", "Chimp", "Bonobo", "Gorilla", "Orangutan"], randomlen=True)
```


---

## 3. Querying tree structure

### Find the MRCA of a set of leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
mrca = tree.find_mrca({"A", "B"})
print(mrca)
```

### Get all descendant leaves below a node

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,((C:1,D:1):1,E:2):1);")
mrca = tree.find_mrca({"C", "E"})
print(tree.remote_children(mrca))
```

### Find leaves near a focal leaf by patristic distance

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
neighbors = tree.nearleafs("A", maxdist=2.1)
print(neighbors)
```

### Find the nearest N leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
print(tree.nearest_n_leaves("A", 3))
```

### Compute pairwise distance between two leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
print(tree.nodedist("A", "D"))
```

### Get the path between two nodes

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
path = tree.nodepath("A", "D")
print(path)
```

### Find diameter and the two farthest leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
maxdist, leaf1, leaf2 = tree.diameter(return_leaves=True)
print(maxdist, leaf1, leaf2)
```

### Find a central or common leaf

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:4);")
leafset = ["A", "B", "C", "D"]
print(tree.find_central_leaf(leafset))
print(tree.find_common_leaf(leafset))
```

---

## 4. Rooting and rerooting

### Midpoint rooting

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5, randomlen=True)
print(tree)
tree.rootmid()
print(tree)
```

### Minimum-variance rooting

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5, randomlen=True)
tree.rootminvar()
```

### Rooting on an outgroup with one taxon

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(Orangutan:4,(Gorilla:3,(Human:2,(Chimpanzee:1,Bonobo:1):1):1):1);")
tree.rootout("Orangutan")
```

### Rooting on an outgroup with multiple taxa

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(Orangutan:4,(Gorilla:3,(Human:2,(Chimpanzee:1,Bonobo:1):1):1):1);")
og = set(["Chimpanzee","Bonobo"])
tree.rootout(og)
```

### Reroot on a specified branch

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
node = tree.find_mrca({"C", "D"})
parent = tree.parent(node)
tree.reroot(parent, node, node1dist=0.4)
```

---

## 5. Modifying topology

### Remove one leaf

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
tree.remove_leaf("D")
```

### Remove several leaves

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);")
tree.remove_leaves(["B", "F"])
```

### Extract a subtree

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
basenode = tree.find_mrca({"A", "B", "C", "D"})
sub = tree.subtree(basenode)
print(sub.newick())
```

### Prune a subtree out of the tree

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
print("Original tree:", tree.newick())
basenode = tree.find_mrca({"E", "F"})
pruned = tree.prune_subtree(basenode)
print("Remaining tree:", tree.newick())
print("Pruned subtree:", pruned.newick())
```

### Insert an internal node

```python
import phylotreelib as pt
from phylotreelib import Branchstruct
tree = pt.Tree.from_string("(A:1,B:1,C:1,D:1);")
newnode = tree.insert_node(0, ["A", "B"], Branchstruct(length=0.5, label="newsplit"))
print(newnode)
print(tree.newick())
```

### Graft one tree onto another

```python
import phylotreelib as pt
left = pt.Tree.from_string("((A:1,B:1):1,C:1);")
right = pt.Tree.from_string("(X:1,Y:1);")
node = left.find_mrca({"A", "B"})
left.graft(right, node1=node, blen1=0.2, blen2=0.3, graftlabel="g2")
print(left.newick())
```

### SPR with specified prune and regraft points

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
prune_node = tree.find_mrca({"C", "D"})
regraft_node = tree.find_mrca({"E", "F"})
tree.spr(prune_node=prune_node, regraft_node=regraft_node)
print(tree.newick())
```

### SPR with specified prune point and random regraft point

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5, randomlen=True)
prune_node = next(iter(tree.leaves))
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

## 6. Representative pruning

### Keep leaves that maximize retained tree length

```python
import phylotreelib as pt
with pt.Nexustreefile("big.tree") as tf:
    tree = tf.readtree()
tree.prune_maxlen(nkeep=50)
```

### Keep leaves that maximize retained tree length while forcing some leaves to remain

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=200, randomlen=True)
tree.prune_maxlen(nkeep=30, keeplist=["s001", "s120", "s178"])
```

### Ask which leaves would be kept without actually pruning, color those on original tree

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=50, randomlen=True)
tree.rootminvar()
keepers = tree.prune_maxlen(nkeep=10, return_leaves=True)
with open("colorprune.nexus", "w") as outfile:
    outfile.write(tree.nexus(colorlist=keepers, colorfg="#FF0000"))
```

---

## 7. Comparing trees

### Robinson–Foulds distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=30, randomlen=True)
t2 = t1.copy_treeobject()
t2.spr()
print(t1.treedist_RF(t2))
print(t1.treedist_RF(t2, normalise=True))
```

### Rooted RF distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=20, randomlen=True)
t2 = t1.copy_treeobject()
t2.rootmid()
print(t1.treedist_RF(t2, rooted=False))
print(t1.treedist_RF(t2, rooted=True))
```

### Path-difference distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=30, randomlen=True)
t2 = t1.copy_treeobject()
for _ in range(3):
    t2.spr()
print(t1.treedist_pathdiff(t2))
```

### Request detailed RF output

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

## 8. Distance matrices

### Build neighbor-joining and UPGMA trees from a nested distance dictionary

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

### Build neighbor-joining tree from a NumPy array

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

## 9. Consensus trees and other tree summaries

**Note**: If you just want a summary tree (consensus/MCC/MBC/HIPSTR) from a set of trees, the command-line tool sumt is usually the simplest option. The examples below show how to do the same thing programmatically with phylotreelib.

* `sumt` GitHub page: [https://github.com/agormp/sumt](https://github.com/agormp/sumt/)

More detailed discussion of summary tree methods is in `docs/api.md`.

### Build a majority-rule consensus tree

```python
import phylotreelib as pt

ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

contree = pt.build_sumtree(ts, treetype="con", blen="biplen", rooting=None)
print(contree.newick())
```

### Build a majority-rule consensus tree with all compatible bipartitions added

```python
import phylotreelib as pt

ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

contree = pt.build_sumtree(ts, treetype="all", blen="biplen", rooting=None)
print(contree.newick())
```

### Build an MCC tree with "common ancestor depths"

```python
import phylotreelib as pt

# Need to keep track of clades and observed topologies to construct MCC
ts = pt.TreeSummary(trackclades=True, tracktopo=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

mcctree = pt.build_sumtree(ts, treetype="mcc", blen="cadepth")
print(mcctree.newick())
```

### Build an MBC tree and mid-point root it

```python
import phylotreelib as pt

# Need to keep track of bipartitions + observed topologies to construct MBC
ts = pt.TreeSummary(trackbips=True, tracktopo=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

mbctree = pt.build_sumtree(ts, treetype="mbc", blen="input", rooting="mid")
print(mbctree.newick())
```


---

## 10. Exporting with labels, metadata, and color

### Plain Newick output

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
with open("tree.newick", "w") as out:
    out.write(tree.newick())
```

### NEXUS output with translate block

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

### Print branch metadata in BEAST/FigTree-style comments

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

### Color leaves selected by a clade query

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

### Color leaves selected by name pattern

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

## 11. Error handling

```python
import phylotreelib as pt
try:
    with pt.Newicktreefile("broken_treefile.newick") as tf:
        tree = tf.readtree()
except pt.TreeError as err:
    print("Tree-related error:", err)
```

---


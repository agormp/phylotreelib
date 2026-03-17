# phylotreelib

[![PyPI version](https://img.shields.io/pypi/v/phylotreelib)](https://pypi.org/project/phylotreelib/)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/phylotreelib?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/phylotreelib)

`phylotreelib` is a Python library for reading, writing, analyzing, manipulating, comparing, and summarizing phylogenetic trees, and is designed to be used in scripts or interactively in notebooks.

![](https://github.com/agormp/phylotreelib/raw/main/treefig.png?raw=true)

Main capabilities include:

- reading and writing Newick and NEXUS tree files
- querying tree structure (MRCA, descendants, paths, distances, neighborhoods)
- rerooting by midpoint, outgroup, or minimum-variance criteria
- pruning leaves or subtrees, regrafting, and other topology edits
- computing tree distances such as Robinson–Foulds and path-difference distance
- building trees from distance matrices (neighbor joining and UPGMA)
- summarizing posterior tree samples into consensus, MCC, MBC, and HIPSTR-style summary trees
- exporting trees with labels, metadata, and optional coloring for visualization in tree viewers (e.g., FigTree)
- library has been optimized for high speed and low memory consumption
- NOTE: labels are interpreted as belonging to branches (bipartitions), not to internal nodes, and this association is maintained after re-rooting etc.

## Installation

```bash
python3 -m pip install phylotreelib
```

Upgrade to the latest version:

```bash
python3 -m pip install --upgrade phylotreelib
```

## Documentation

This repository uses three main documentation files:

- `README.md` — (this file) overview and quick start
- `docs/recipes.md` — copy/paste-ready examples and workflows
- `docs/api.md` — more extensive reference to the API, plus advanced notes

The README gives examples of the most common workflows. More specialized details, and the theory behind some of the summary-tree machinery, are in the `docs/` files.

## Terminology note: “depth” vs “height” (tree orientation)

In **phylotreelib**, trees are mostly treated as *rooted at the bottom and having leaves at the top*.
Accordingly, it uses:

- **depth(node)** = distance from the **tips (leaves)** to the node (i.e., “time before most recent leaf”)

Similarly terms like above/below and upper/lower are based on the same perceived orientation of the tree.

However, for reasons clear only to computer scientists, in many programming libraries and APIs, trees are often treated as rooted at the top, and the same concept is referred to as **height**. I am in the process of changing the phylotreelib names and documentation to follow that standard practice. The deprecated depth-names will stay functional for one major release cycle after changing them (in the upcoming phylotreelib version 2.0).


---

## Quick start

### Construct a tree from a Newick string, print tabular representation

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(Orangutan:4,(Gorilla:3,(Human:2,(Chimpanzee:1,Bonobo:1):1):1):1);")
print(tree)
```

Output:

```text
+---------------------------------------+
|  parent  |    child     |  branchlen  |
+---------------------------------------+
|       0  |           1  |          1  |
|       0  |   Orangutan  |          4  |
|       1  |           2  |          1  |
|       1  |     Gorilla  |          3  |
|       2  |           3  |          1  |
|       2  |       Human  |          2  |
|       3  |      Bonobo  |          1  |
|       3  |  Chimpanzee  |          1  |
+---------------------------------------+

5 Leaves:
----------
Bonobo
Chimpanzee
Gorilla
Human
Orangutan
```

### Read one tree from file, reroot it, print all root-to-tip distances

```python
import phylotreelib as pt
with pt.Nexustreefile("mytreefile.nexus") as treefile:
    tree = treefile.readtree()
tree.rootminvar()
for tip in sorted(tree.leaves):
    dist = tree.nodedist(tree.root, tip)
    print(f"{tip:<12s} {dist:.3f}")
```

### Query tree structure

The examples below illustrate some of the most useful tree queries: finding an MRCA, listing descendants, computing patristic distance, finding nearby leaves by patristic distance.

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(Orangutan:4,(Gorilla:3,(Human:2,(Chimpanzee:1,Bonobo:1):1):1):1);")
mrca = tree.find_mrca({"Human", "Chimpanzee"})
print("MRCA:", mrca)
print("Descendant leaves:", tree.remote_children(mrca))
print("Distance Human–Bonobo:", tree.nodedist("Human", "Bonobo"))
print("Leaves within distance 4.5 of Human:", tree.nearleafs("Human", 4.5))
print("Nearest 3 leaves to Human:", tree.nearest_n_leaves("Human", 3))
```

### Remove selected leaves, keep rest of topology

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,((C:1,D:1):1,(E:1,F:1):1):1);")
tree.remove_leaves(["B", "D"])
print(tree.newick())
```

Output:

```text
(((E:1,F:1):1,C:2):1,A:2);
```

### Perform a specified SPR move

```python
import phylotreelib as pt
tree = pt.Tree.from_string("(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);")
prune_node = tree.find_mrca({"C", "D"})
regraft_node = tree.find_mrca({"E", "F"})
tree.spr(prune_node=prune_node, regraft_node=regraft_node)
print(tree.newick())
```

### Perform a semi-random or fully random SPR move

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=12, randomlen=True)
# semi-random: choose prune node, randomize regraft point
prune_node = next(iter(tree.intnodes))
tree.spr(prune_node=prune_node)
# fully random
for _ in range(3):
    tree.spr()
```

### Summarize posterior tree samples

If you are working with posterior tree samples (e.g. from BEAST/MrBayes), you can build a summary tree in a few explicit steps:

1) **TreeSummary** accumulates frequencies and (optionally) branch-length/node-depth statistics
2) **SummaryTreeBuilder** constructs a summary-tree topology
3) **TreePostProcessor** adds attributes to the tree annotating it with e.g. support values / branch-length stats
4) **PrintSpec** controls what attributes are written as branch labels and/or Nexus meta-comments

For command-line use, the separate tool [`sumt`](https://github.com/agormp/sumt) wraps these steps in a single command.

```python
import phylotreelib as pt

# 1) Accumulate bipartition frequencies + mean branch lengths
ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

# 2) Build a majority-rule consensus topology
stb = pt.SummaryTreeBuilder(ts)
sumtree = stb.contree(cutoff=0.5, allcompat=False)

# 3) Annotate/paint the tree from the collected summaries
tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(sumtree)

# 4) Configure printing (labels + optional meta-comments)
pt.configure_sumtree_printing(sumtree, treetype="con", blen="biplen", print_meta=True, precision=7)

print(sumtree.newick())
```

### Compute and export a summary tree in one step

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
    trackci=True,
)

with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)
sumtree = pt.build_sumtree(ts, treetype="mcc", blen="meandepth", rooting=None,)

pt.configure_sumtree_printing(
	sumtree,
	treetype="mcc",
	blen="meandepth",
	trackci=True,
	print_meta=True
)

with open("posterior_summary.nexus", "w") as out:
    out.write(sumtree.nexus())
```

### Save in different formats and include metadata

The printing helpers let you choose labels, branch lengths, and metadata comments independently of the tree itself.

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
# attach a branch attribute that can later be printed
mrca = tree.find_mrca({"A", "B"})
parent = tree.parent(mrca)
tree.set_branch_attribute(parent, mrca, "posterior", 0.97)
pt.configure_basic_printing(
    tree,
    print_meta=True,
    branch_attrs=("posterior",),
    precision=6,
)
print(tree.newick())
print(tree.nexus(translateblock=True))
```

### Color a selected set of leaves in NEXUS output

This is useful for visualization in for instance FigTree. Here, the leaves to color are obtained by a biological query on the tree (specifically: the descendants of the MRCA of the mink leaves).

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((Human_A:1,Human_B:1):1,(Mink_A:1,Mink_B:1):1,(Bat_A:1,Bat_B:1):1);")
mink_mrca = tree.find_mrca({"Mink_A", "Mink_B"})
color_me = sorted(tree.remote_children(mink_mrca))
with open("colored.nexus", "w") as outfile:
    outfile.write(tree.nexus(translateblock=False, colorlist=color_me,
							 colorfg="#FF0000", colorbg="#000000"))
```


---

## Core concepts

### `Tree` class

`Tree` is the central class. Most workflows involve reading or constructing one or more `Tree` objects and then querying or modifying them.

Useful attributes:

- `tree.leaves` — set of leaf names
- `tree.intnodes` — set of internal node IDs
- `tree.nodes` — set of all nodes
- `tree.root` — current root node ID

A tabular representation can be obtained with:

```python
print(tree)
```

### Tree files

Recommended usage:

```python
import phylotreelib as pt
with pt.Treefile("trees.nex") as tf:
    tree = tf.readtree()
```

Available classes:

- `Treefile` — file-format autodetection
- `Newicktreefile`
- `Nexustreefile`

### Tree sets and tree summaries

- `TreeSet` stores multiple trees with the same leaves.
- `TreeSummary` accumulates summary statistics across trees and constructs summary trees such as consensus, MCC, MBC, and HIPSTR.

---

## Where to go next

- For more worked examples, see [`docs/recipes.md`](docs/recipes.md)
- For the API reference and advanced notes, see [`docs/api.md`](docs/api.md)


### Related tools

- [`sumt`](https://github.com/agormp/sumt) — command-line interface for summary-tree workflows built on `phylotreelib`

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


## Version 2.0.0

Version 2 has recently been released, and contains a number of changes to the API compared to version 1.x.x.

- See the **Quick start** below for updated examples.
- For details on API changes and how to update existing scripts, see **[Upgrading from 1.x to 2.x](#upgrading-from-1x-to-2x)**.
- For a feature overview, see **[Changelog (2.0.0)](#changelog-200)**.


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


---

## Quick start

### Terminology note: “depth” vs “height” (tree orientation)

In **phylotreelib**, trees are mostly treated as *rooted at the bottom and having leaves at the top*.
Accordingly, it uses:

- **depth(node)** = distance from the **tips (leaves)** to the node (i.e., “time before most recent leaf”)

Similarly terms like above/below and upper/lower are based on the same perceived orientation of the tree.

However, for reasons clear only to computer scientists, in many programming libraries and APIs, trees are often treated as rooted at the top, and the same concept is referred to as **height**. I am in the process of changing the phylotreelib names and documentation to follow that standard practice. The deprecated depth-names will stay functional for one major release cycle after changing them (for the duration of version 2).


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


---

## Changelog (2.0.0)

### Added

- **Explicit summary-tree pipeline** via `SummaryTreeBuilder` and `TreePostProcessor` (with optional `CADepthEstimator` for CA-depths and credible intervals).
- **Credible interval support** via `QuantileAccumulator` and `trackci` / `ci_probs` in `TreeSummary`.
- **PrintSpec** plus helpers `configure_basic_printing(...)` and `configure_sumtree_printing(...)` for consistent Newick/NEXUS output.

### Changed

- `Tree.graft(...)` now grafts onto a specified edge `(parent -> child)` and uses `connect_length` / `connect_branchstruct` for the connecting branch.
- `Tree.spr(...)` updated to match the new grafting model and (by default) preserve total tree length.
- `Tree.subtree(...)` / `Tree.prune_subtree(...)` now return `(subtree, basal_branch)`; leaf subtrees return the leaf name as `str`.
- `Tree.reroot(...)` signature updated to use the same distance/fraction placement convention as other editing methods.
- `Tree.newick(...)` / `Tree.nexus(...)` now treat `None` as “use PrintSpec defaults”.

### Deprecated

- `Tree.deroot()` → use `Tree.collapse_bifurcating_root()`.

### Fixed / improved

- Multiple refactorings leading to performance improvements (speed and memory) and bug fixes.


---

## Upgrading from 1.x to 2.x

phylotreelib 2.0 introduces a few intentional API breaks to make tree editing, summary trees, and printing more explicit and less “magic”.

If you maintain scripts that call phylotreelib directly, the sections below show the key signature/return changes and how to update call sites.

### Summary trees

`TreeSummary.compute_sumtree(...)` is still available as a convenience wrapper (so many 1.x workflows continue to work), but the preferred style is now the explicit pipeline:

- `TreeSummary` --> `SummaryTreeBuilder` --> `TreePostProcessor` (+ optional `CADepthEstimator`) --> `PrintSpec`

See the “Summarize posterior tree samples” quick-start example above, and the worked examples in `docs/recipes.md`.

### Tree.graft: signature change

**Old (1.x):**

```python
Tree.graft(self, other, node1, node2=None, blen1=0, blen2=0,
           graftlabel=None, graft_with_other_root=False)
```

Meaning (1.x): attach `other` below `node1` (possibly re-rooting `other` around `node2`), by inserting a new internal node above `node1` with length `blen1`, and connecting `other` above that with length `blen2` (or optionally skip the extra basal branch).

**New (2.x):**

```python
Tree.graft(self, other, parent, child,
           dist_from_parent=None, frac_from_parent=None,
           connect_length=0.0, connect_branchstruct=None, graftlabel=None)
```

Meaning (2.x): insert the graft point on an *existing edge* `(parent -> child)` and then attach `other` (or a single leaf name) below that graft point via a *separate connecting branch*.

**How to convert typical 1.x calls**

1) **Grafting “below node1” (typical old usage)**

Old:

```python
tree.graft(other, node1=X, blen1=blen1, blen2=blen2)
```

New (equivalent intent):

```python
parent = tree.parent(X)
tree.graft(
    other,
    parent=parent, child=X,
    # place graftpoint on the incoming edge to X:
    dist_from_parent=blen1,           # if you want a specific distance
    # or frac_from_parent=0.5,        # if you just want “midway”
    connect_length=blen2,
)
```

2) **Old `graft_with_other_root=True`**

In 2.x, the caller should explicitly decide what `other` is before calling `graft(...)`:
- If you want to attach a *single leaf*, pass a string `other="TaxonName"` (and skip the “extra basal branch” complexity entirely).
- If you want to attach a Tree with a specific root, re-root `other` first (e.g. `other.reroot(...)`) and then call `graft(...)`.

### Tree.spr: signature change (and length preservation)

**Old (1.x):**

```python
Tree.spr(self, prune_node=None, regraft_node=None)
```

**New (2.x):**

```python
Tree.spr(self, prune_node=None, regraft_node=None,
         dist_from_parent=None, frac_from_parent=None)
```

New behaviour highlights:

- Regrafting happens by inserting a graft point on the incoming edge to `regraft_node`.
- If `dist_from_parent` and `frac_from_parent` are both `None`, the placement is randomized (uniform frac in (0,1)).
- Total tree length is preserved by default (the regraft edge is split; the pruned basal branch length is reused as the connecting branch length).

**How to convert**

Old:

```python
tree.spr(prune_node=P, regraft_node=R)
```

New (same defaults):

```python
tree.spr(prune_node=P, regraft_node=R)
```

New (explicit placement on incoming edge to `regraft_node`):

```python
parent = tree.parent(R)
tree.spr(prune_node=P, regraft_node=R, frac_from_parent=0.25)
# (places the graftpoint 25% down the edge parent->R)
```

### Tree.subtree / Tree.prune_subtree: return value change

**Old (1.x):**

```python
Tree.subtree(self, basenode, return_basalbranch=False) -> Tree
Tree.prune_subtree(self, basenode) -> Tree
```

**New (2.x):**

```python
Tree.subtree(self, basenode) -> (subtree, basal_branch)
Tree.prune_subtree(self, basenode) -> (subtree, basal_branch)
```

Notes (2.x):

- `basal_branch` is a **copy** of the Branchstruct on the incoming edge to `basenode`.
- If `basenode` is a leaf: `subtree` is the **leaf name** (`str`), not a mini tree.

**How to convert**

Old:

```python
sub = tree.subtree(node)
```

New:

```python
sub, basal = tree.subtree(node)
```

Old:

```python
sub = tree.prune_subtree(node)
```

New:

```python
sub, basal = tree.prune_subtree(node)
# basal.length is the length removed when detaching the subtree
```

### Tree.reroot: signature change

**Old (1.x):**

```python
Tree.reroot(self, node1, node2=None, polytomy=False, node1dist=0.0)
```

**New (2.x):**

```python
Tree.reroot(self, node1, node2=None, polytomy=False,
            dist_from_node1=None, frac_from_node1=0.5)
```

Conversion:

- If you previously used `node1dist`, rename it to `dist_from_node1`.
- If you previously relied on the default placement, you can now do that via `frac_from_node1` (default 0.5).

### Root-collapsing helper rename

**Old (1.x):**

```python
Tree.deroot(self)
```

**New (2.x):**

```python
Tree.collapse_bifurcating_root(self)
```

(And `Tree.deroot()` is kept as a deprecated alias.)

### Printing: PrintSpec defaults

**Old (1.x):**

```python
Tree.newick(self, printdist=True, printlabels=True, labelfield='label', precision=6, ...)
Tree.nexus(self, printdist=True, printlabels=True, labelfield='label', precision=6, ...)
```

**New (2.x):**

```python
Tree.newick(self, printdist=None, printlabels=None, labelfield=None, precision=None, ..., print_meta=None)
Tree.nexus(self,  printdist=None, printlabels=None, labelfield=None, precision=None, ..., print_meta=None)
```

In 2.x, passing `None` means “use the Tree’s current PrintSpec default”. This makes it easy to configure printing once:

```python
pt.configure_basic_printing(tree, precision=7, print_meta=True, branch_attrs=("posterior",))
print(tree.newick())          # uses PrintSpec
print(tree.nexus())           # uses PrintSpec
```

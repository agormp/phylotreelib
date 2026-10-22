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


## Version 2

Version 2 introduced several API changes compared with version 1.x.

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

This repository uses four main documentation files:

- `README.md` — (this file) overview and quick start
- `docs/recipes.md` — copy/paste-ready examples and workflows
- `docs/api.md` — more extensive reference to the API, plus advanced notes
- `CHANGELOG.md` — version history and migration notes

The README gives examples of the most common workflows. More specialized details, and the theory behind some of the summary-tree machinery, are in the `docs/` files.


---

## Quick start

### Terminology note: tree orientation and "height" vs "depth"

**Note for users upgrading from version 1:** In **phylotreelib** v. 2, trees are now treated as *having the root at the top and leaves at the bottom* — matching the standard convention in most phylogeny software and computer science.

Accordingly, it uses:

- **height(node)** = distance from the **tips (leaves)** to the node (i.e., "time before most recent leaf")

In v1.x, this concept was called **depth**. All depth-based names were renamed to height-based names in v2.0.0. See [Upgrading from 1.x to 2.x](#upgrading-from-1x-to-2x) for the full rename table.

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
Tree length: 14.0

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
print(tree.newick(printdist=False))
# semi-random: choose prune node, randomize regraft point
prune_node = next(iter(tree.possible_spr_prune_nodes()))
tree.spr(prune_node=prune_node)
print(tree.newick(printdist=False))
# fully random
for _ in range(10):
    tree.spr()
    print(tree.newick(printdist=False))
```

### Summarize posterior tree samples

If you are working with posterior tree samples (e.g. from BEAST/MrBayes), you can build a summary tree in a few explicit steps:

1) **TreeSummary** accumulates frequencies and (optionally) branch-length/node-height statistics
2) **SummaryTreeBuilder** constructs a summary-tree topology
3) **TreePostProcessor** adds attributes to the tree annotating it with e.g. support values / branch-length stats
4) **Tree.newick() / Tree.nexus()** arguments select branch labels and Nexus meta-comments for each output call

For command-line use, the separate tool [`sumt`](https://github.com/agormp/sumt) wraps these steps in a single command.

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/data/primate-mtDNA.trees"
urllib.request.urlretrieve(URL, "primate-mtDNA.trees")

import phylotreelib as pt

# 1) Accumulate bipartition frequencies + mean branch lengths
ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

# 2) Build a majority-rule consensus topology
stb = pt.SummaryTreeBuilder(ts)
sumtree = stb.contree(cutoff=0.5, allcompat=False)

# 3) Annotate/paint the tree from the collected summaries
tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(sumtree)

# 4) Select labels and metadata in the output call
print(sumtree.newick(
    labelfield="bipartition_cred",
    branch_attributes=("bipartition_cred", "length", "length_sd"),
    print_meta=True,
    precision=7,
))
```

### Compute credible intervals and export a summary tree

This small, self-contained example computes an MCC tree with mean clade heights and 80% and 95% credible intervals. For a real analysis, feed trees from `Nexustreefile` into the same accumulator, excluding burn-in first.

```python
import phylotreelib as pt

samples = (
    "((A:1,B:1):2,(C:2,D:2):1);",
    "((A:1.5,B:1.5):2.5,(C:2,D:2):2);",
    "((A:1.2,B:1.2):2.8,(C:2.5,D:2.5):1.5);",
)
ts = pt.TreeSummary(
    trackclades=True,
    trackheight=True,
    tracktopo=True,
    trackci=True,
    ci_probs=(0.8, 0.95),
)
for sample in samples:
    ts.add_tree(pt.Tree.from_string(sample))

sumtree = pt.build_sumtree(ts, treetype="mcc", blen="cladeheight")

# All populated annotations, including the intervals, are exported automatically.
with open("posterior_summary.nexus", "w") as out:
    out.write(sumtree.nexus(print_meta=True))

# A smaller export: heights, medians, and only the 95% interval.
print(sumtree.newick(
    print_meta=True,
    node_attributes=("height", "height_median", "ci"),
    branch_attributes=(),
    ci_labels=("95%_CI",),
))
```

Specify both `trackci=True` and `ci_probs`: the former enables quantile tracking, while the latter requests interval probabilities. These are central, equal-tailed credible intervals, not HPD intervals. The Nexus metadata contains `height`, `height_median`, `height_80%_CI={lower,upper}`, and `height_95%_CI={lower,upper}`, including annotations on the root. Variable-height quantiles are approximate; effectively constant node heights retain the observed value for the median and both interval endpoints.

See [the recipes](docs/recipes.md#9-consensus-trees-and-other-tree-summaries) for common-ancestor-height intervals and [the API notes](docs/api.md#106-what-the-tracked-quantiles-are-used-for) for branch-length intervals.

### Save in different formats and include metadata

Pass printing options directly to `newick()` or `nexus()`. Options apply only to that call. An ordinary dictionary can be reused for several calls:

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
# attach a branch attribute that can later be printed
mrca = tree.find_mrca({"A", "B"})
parent = tree.parent(mrca)
tree.set_branch_attribute(parent, mrca, "posterior", 0.97)
print_options = {
    "print_meta": True,
    "node_attributes": (),
    "branch_attributes": ("posterior",),
    "precision": 6,
}
print(tree.newick(**print_options))
print(tree.nexus(translateblock=True, **print_options))
print(tree.newick())  # Plain output; previous calls do not change its defaults
```

### Color a selected set of leaves in NEXUS output

This is useful for visualization in for instance FigTree. Here, the leaves to color are obtained by a biological query on the tree (specifically: the descendants of the MRCA of the mink leaves).

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((Human_A:1,Human_B:1):1,(Mink_A:1,Mink_B:1):1,(Bat_A:1,Bat_B:1):1);")
mink_mrca = tree.find_mrca({"Mink_A", "Mink_B"})
color_me = sorted(tree.remote_children(mink_mrca))
print(tree.nexus(translateblock=False, colorlist=color_me,
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
- `TreeSummary` accumulates summary statistics across trees.
- `SummaryTreeBuilder` constructs summary trees such as consensus, MCC, MBC, and HIPSTR.
- `TreePostProcessor` adds attributes to the tree annotating it with e.g. support values / branch-length stats.

---

## Where to go next

- For more worked examples, see [`docs/recipes.md`](docs/recipes.md)
- For the API reference and advanced notes, see [`docs/api.md`](docs/api.md)


### Related tools

- [`sumt`](https://github.com/agormp/sumt) — command-line interface for summary-tree workflows built on `phylotreelib`


---

## Changelog (2.0.0)

These historical notes describe version 2.0.0. The unreleased printing API replaces its stored printing settings; see [Printing without stored defaults](#printing-without-stored-defaults-unreleased).

### Added

- **Explicit summary-tree pipeline** via `SummaryTreeBuilder` and `TreePostProcessor` (with optional `CAHeightEstimator` for CA-heights and credible intervals).
- **Credible interval support** via `QuantileAccumulator` and `trackci` / `ci_probs` in `TreeSummary`.
- **PrintSpec** plus helpers `configure_basic_printing(...)` and `configure_sumtree_printing(...)` for consistent Newick/NEXUS output.

### Changed

- `Tree.graft(...)` now grafts onto a specified edge `(parent -> child)` and uses `connect_length` / `connect_branchstruct` for the connecting branch.
- `Tree.spr(...)` updated to match the new grafting model and (by default) preserve total tree length.
- `Tree.subtree(...)` / `Tree.prune_subtree(...)` now return `(subtree, basal_branch)`; leaf subtrees return the leaf name as `str`.
- `Tree.reroot(...)` signature updated to use the same distance/fraction placement convention as other editing methods.
- `Tree.newick(...)` / `Tree.nexus(...)` now treat `None` as "use PrintSpec defaults".

### Renamed (terminology: `depth` → `height`)

All names relating to node heights — the distance from the tips to an internal node, i.e. "time before the most recent leaf" — have been renamed from `depth`-based to `height`-based. This aligns with the convention used in most phylogenetic software and programming libraries. **This is a breaking change; the old names are no longer available.**

| v1.x name | v2.x name |
|---|---|
| `TreeSummary(..., trackdepth=True, ...)` | `TreeSummary(..., trackheight=True, ...)` |
| `Tree.nodedepth(node)` | `Tree.nodeheight(node)` |
| `Tree.set_blens_from_depths()` | `Tree.set_blens_from_heights()` |
| `TreePostProcessor.set_mean_node_depths(tree)` | `TreePostProcessor.set_mean_node_heights(tree)` |
| `CADepthEstimator` | `CAHeightEstimator` |
| `build_sumtree(..., blen="meandepth", ...)` | `build_sumtree(..., blen="cladeheight", ...)` |
| `build_sumtree(..., blen="cadepth", ...)` | `build_sumtree(..., blen="caheight", ...)` |
| Node attribute `depth` | Node attribute `height` |
| Node attribute `depth_sd` | Node attribute `height_sd` |
| Node attribute `depth_median` | Node attribute `height_median` |

### Deprecated

- `Tree.deroot()` → use `Tree.collapse_bifurcating_root()`.

### Fixed / improved

- Multiple refactorings leading to performance improvements (speed and memory) and bug fixes.


---

## Upgrading from 1.x to 2.x

phylotreelib 2.0 introduces a few intentional API breaks to make tree editing, summary trees, and printing more explicit and less "magic".

If you maintain scripts that call phylotreelib directly, the sections below show the key signature/return changes and how to update call sites.

### Terminology rename: `depth` → `height`

In v1.x, the distance from the tips to an internal node was called **depth**. In v2.0.0 this has been renamed to **height** throughout the library. The rename is purely mechanical: a search-and-replace across your scripts using the table in the [Changelog](#changelog-200) above is sufficient.

The affected call sites are:

- `TreeSummary` constructor: `trackdepth=` → `trackheight=`
- `build_sumtree` and `configure_sumtree_printing`: `blen="cladedepth"` → `"cladeheight"`, `blen="cadepth"` → `"caheight"`
- `Tree` methods: `nodedepth()` → `nodeheight()`, `set_blens_from_depths()` → `set_blens_from_heights()`
- Node attributes on summary trees: `depth` / `depth_sd` / `depth_median` → `height` / `height_sd` / `height_median`

### Summary trees

`TreeSummary.compute_sumtree(...)` is still available as a convenience wrapper (so many 1.x workflows continue to work), but the preferred style is now the explicit pipeline:

- `TreeSummary` --> `SummaryTreeBuilder` --> `TreePostProcessor` (+ optional `CAHeightEstimator`) --> `Tree.newick()` / `Tree.nexus()`

See the "Summarize posterior tree samples" quick-start example above, and the worked examples in `docs/recipes.md`.

### Tree.graft: signature change

**Old (1.x):**

```python
Tree.graft(self, other, node1, node2=None, blen1=0, blen2=0,
           graftlabel=None, graft_with_other_root=False)
```

Meaning (1.x): attach `other` on incoming branch to `node1` (possibly re-rooting `other` around `node2`): insert a new internal node on the incoming branch to `node1` with length `blen1`, and connect `other` below that (meaning closer to the tips) with length `blen2` (or optionally skip the extra basal branch).

**New (2.x):**

```python
Tree.graft(self, other, parent, child,
           dist_from_parent=None, frac_from_parent=None,
           connect_length=0.0, connect_branchstruct=None, graftlabel=None)
```

Meaning (2.x): insert the graft point on an *existing edge* `(parent -> child)` and then attach `other` (or a single leaf name) below that graft point via a *separate connecting branch*.

**How to convert typical 1.x calls**

1) **Grafting on incoming edge to node1 (typical old usage)**

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
    # or frac_from_parent=0.5,        # if you just want "midway"
    connect_length=blen2,
)
```

2) **Old `graft_with_other_root=True`**

In 2.x, the caller should explicitly decide what `other` is before calling `graft(...)`:
- If you want to attach a *single leaf*, pass a string `other="TaxonName"`.
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
- Total tree length is preserved by default (the regraft edge is split according to `dist_from_parent` or `frac_from_parent`; the length of the incoming basal branch to the prune_node is reused as the connecting branch length).

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
# parent = tree.parent(R)
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

### Printing without stored defaults (unreleased)

To simplify the workflow, printing is again configured where the tree is exported: in `newick()` or `nexus()`. The arguments at that call describe the output; there is no earlier configuration step or stored formatting state to keep track of. An ordinary options dictionary can still be reused across calls.

This is a breaking API change after 2.3.1. `PrintSpec`, `configure_basic_printing()`, `configure_sumtree_printing()`, `Tree.set_print_spec()`, and `Tree.get_print_spec()` have been removed. Pass the options directly to each output call.

**Before (2.0.0 through 2.3.1):**

```python
pt.configure_basic_printing(tree, precision=7, print_meta=True, branch_attrs=("posterior",))
print(tree.newick())
print(tree.nexus())
```

**Now:**

```python
print(tree.newick(precision=7, print_meta=True, node_attributes=(), branch_attributes=("posterior",)))
print(tree.nexus(precision=7, print_meta=True, node_attributes=(), branch_attributes=("posterior",)))
```

The helper arguments `node_attrs` and `branch_attrs` become `node_attributes` and `branch_attributes` on the output methods. For summary trees, `sumtree.nexus(print_meta=True)` is enough to include all populated annotations; knowing their individual names is optional.

`Tree.newick()` and `Tree.nexus()` default to printing lengths and labels, using `labelfield="label"`, `precision=7`, and `print_meta=False`. `TreeSet` uses the same defaults, including seven significant digits (previously six).

With `print_meta=True`, attribute selections work independently for nodes and branches:

| Selection | Output |
|---|---|
| Omitted or `None` | All populated annotation fields, including stored intervals |
| `()` or `[]` | No metadata of that kind |
| A sequence of names | Only those attributes; include `"ci"` to emit intervals |

Internal accumulators and temporary branch-height caches are excluded from automatic output. Missing/`None` values and empty strings or collections are skipped; zero and `False` are retained. `ci_labels=None` includes all stored intervals, an empty sequence suppresses intervals, and a sequence of labels selects those present.

An explicit `None` for a scalar option uses its fixed default. `print_meta=False` always suppresses node and branch metadata, even if attributes are selected. These settings never persist on the tree. When migrating selective exports, use empty tuples for unwanted metadata and include `"ci"` in explicit attribute selections that should print intervals.

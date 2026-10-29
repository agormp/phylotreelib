# phylotreelib

[![PyPI version](https://img.shields.io/pypi/v/phylotreelib)](https://pypi.org/project/phylotreelib/)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/phylotreelib?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/phylotreelib)

`phylotreelib` is a Python library for reading, writing, analyzing, comparing and manipulating phylogenetic trees. It can be used in scripts or interactively in notebooks.

![](https://github.com/agormp/phylotreelib/raw/main/treefig.png?raw=true)

Main capabilities include:

- reading and writing Newick and Nexus trees, including tree annotations
- querying tree structure, MRCAs, paths and distances
- rooting, pruning, grafting and other topology edits
- comparing trees using Robinson–Foulds and path-difference distances
- building neighbor-joining and UPGMA trees from distance matrices
- summarizing tree samples as consensus, MCC, MBC and HIPSTR trees

## Installation

This README describes phylotreelib 3.0.1, which requires Python 3.10 or newer.

```bash
python3 -m pip install phylotreelib
```

To upgrade an existing installation:

```bash
python3 -m pip install --upgrade phylotreelib
```

## Documentation

- This README gives a quickstart and common examples.
- [Recipes](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md) give worked examples for reading, editing, comparing and summarizing trees.
- The [API reference](https://github.com/agormp/phylotreelib/blob/main/docs/api.md) describes arguments, data fields and copying/ownership rules.
- The [changelog](https://github.com/agormp/phylotreelib/blob/main/CHANGELOG.md) records notable changes.

If you are updating scripts, see the [3.0 migration guide](https://github.com/agormp/phylotreelib/blob/main/docs/api.md#upgrading-to-30) and the [older 1.x upgrade examples](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#13-upgrading-existing-1x-scripts).

## Quick start

### Construct and query a tree

```python
import phylotreelib as pt

tree = pt.Tree.from_string("(Orangutan:4,(Gorilla:3,(Human:2,(Chimpanzee:1,Bonobo:1):1):1):1);")
print(tree)
print(tree.newick())

mrca = tree.find_mrca({"Human", "Chimpanzee"})
print("Descendant leaves:", tree.remote_children(mrca))
print("Human–Bonobo distance:", tree.nodedist("Human", "Bonobo"))
```

Node height means distance from the most recent tips to the node. Version 1 called this depth. See the [upgrade examples](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#13-upgrading-existing-1x-scripts) for renamed methods and arguments.

### Read trees from a file

`Treefile` detects Newick or Nexus format. Iterate over a file to process trees without storing them all:

```python
import phylotreelib as pt

with pt.Treefile("posterior.trees") as treefile:
    treefile.skip_trees(100)  # Discard burn-in.
    for tree in treefile:
        print(tree.length())
```

Use `Newicktreefile` or `Nexustreefile` when the format is known. `readtree()` returns one Tree, or None at the end of the file. More input examples, including raw tree text, are in the [reading recipes](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#1-reading-trees).

### Root and prune a tree

```python
import phylotreelib as pt

tree = pt.Tree.from_string("((A:1,B:1):1,((C:1,D:1):1,(E:1,F:1):1):1);")
tree.rootmid()
tree.remove_leaves(["B", "D"])
print(tree.newick())
```

See the recipes for [other rooting methods](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#4-rooting-and-rerooting), [grafting](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#graft-one-tree-onto-another) and [SPR moves](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#spr-with-specified-prune-and-regraft-points).

### Read and write annotations

Use `parse_meta=True` to read annotation values. Default reading skips optional annotations but retains supported rootedness declarations and root labels/stem lengths.

```python
import phylotreelib as pt

tree = pt.Tree.from_string('[&R](A[&location="Denmark"]:1,B:1);', parse_meta=True)
print(tree.nodedict["A"].extra_meta["location"])
print(tree.nexus(print_meta=True))
```

Printing options apply to each call. `print_meta=True` includes populated annotations. Select particular keys with `node_attributes` and `branch_attributes`, or use an empty tuple to omit that kind of annotation. Output uses full stored precision unless a numeric `precision` is requested.

`tree.rooted` is True, False or None. Nexus output declares known state with `[&R]` or `[&U]`. Plain Newick does not include a declaration. See the [input and output API](https://github.com/agormp/phylotreelib/blob/main/docs/api.md#2-tree-input-and-output) for positional annotation rules, missing values and output controls.

### Summarize posterior tree samples

`TreeSummary` collects frequencies and optional length/height statistics. `build_sumtree()` constructs and annotates the requested summary tree.

For command-line use, the separate tool [sumt](https://github.com/agormp/sumt) provides these workflows in a single command.

### Compute credible intervals and export a summary tree

This self-contained example builds an MCC tree with mean clade heights and 80% and 95% credible intervals:

```python
import phylotreelib as pt

samples = (
    "((A:1,B:1):2,(C:2,D:2):1);",
    "((A:1.5,B:1.5):2.5,(C:2,D:2):2);",
    "((A:1.2,B:1.2):2.8,(C:2.5,D:2.5):1.5);",
)
summary = pt.TreeSummary(
    trackclades=True,
    trackheight=True,
    tracktopo=True,
    trackci=True,
    ci_probs=(0.8, 0.95),
)
for sample in samples:
    summary.add_tree(pt.Tree.from_string(sample))

sumtree = pt.build_sumtree(summary, treetype="mcc", blen="cladeheight")
print(sumtree.nexus(print_meta=True))
```

For a real analysis, read samples from a file and exclude burn-in before adding observations. Specify both `trackci=True` and `ci_probs`. These are central, equal-tailed intervals. Variable-height quantiles are approximate, while effectively constant node heights retain the observed value.

See the [summary-tree recipes](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#9-consensus-trees-and-other-tree-summaries) for consensus trees, common-ancestor heights and selected interval output. If you merge partial summaries or retain input representatives, follow the [ownership guidance](https://github.com/agormp/phylotreelib/blob/main/docs/api.md#5-summary-ownership).

## More examples

- [Extract and prune subtrees](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#extract-a-subtree)
- [Compare trees](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#7-comparing-trees)
- [Build trees from distance matrices](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#8-distance-matrices)
- [Export labels, annotations and colors](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#10-exporting-with-labels-metadata-and-color)
- [Parsimony calculations](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#12-parsimony)

## Upgrading existing scripts

Version 3.0 removes stored printing settings, reader `strip_comments`, `below_root` and public `.ci` storage. Printing and interval selection are now explicit. Read the [migration table](https://github.com/agormp/phylotreelib/blob/main/docs/api.md#upgrading-to-30) before updating callers.

The [older upgrade recipes](https://github.com/agormp/phylotreelib/blob/main/docs/recipes.md#13-upgrading-existing-1x-scripts) cover the 1.x-to-2.x height renames, grafting/SPR arguments and subtree return values.

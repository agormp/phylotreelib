## phylotreelib recipes

This file contains task-oriented examples for common `phylotreelib` workflows.

The recipes are intentionally practical and somewhat redundant. The goal is to provide blocks of code that can be copied and easily adapted for related tasks.

### Terminology note: tree orientation and "height" vs "depth"

**Note for users upgrading from version 1:** In **phylotreelib** v. 2, trees are now treated as *having the root at the top and leaves at the bottom* — matching the standard convention in most phylogeny software and computer science.

Accordingly, it uses:

- **height(node)** = distance from the **tips (leaves)** to the node (i.e., "time before most recent leaf")

In v1.x, this concept was called **depth**. All depth-based names were renamed to height-based names in v2.0.0. See [older upgrade examples](#13-upgrading-existing-1x-scripts) for the full rename table.

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

#### Import selected positional annotations and retain raw text explicitly

```python
import phylotreelib as pt
source = '[&R](A[&height=1.0,location="north coast"]:0.5[&rate=0.7],B:2)Root:0;'
with pt.Newicktreefile(filecontent=source, parse_meta=True) as tf:
    raw = tf.read_tree_text()  # Complete original description, no Tree construction.
assert raw == source
tree = pt.Tree.from_string(raw, parse_meta=True)
assert tree.rooted is True
assert tree.nodedict["A"].height == 1.0
assert tree.nodedict["A"].extra_meta["location"] == "north coast"
assert tree.child_dict[tree.root]["A"].extra_meta["rate"] == 0.7
assert tree.root_suffix.label == "Root" and tree.root_suffix.length == 0.0
print(tree.nexus(print_meta=True, node_attributes=("height", "location"),
                 branch_attributes=("rate",)))
```

Annotations before a colon belong to the node; those after it belong to its incoming branch, or the root's stem. Default `parse_meta=False` keeps supported declarations and structural suffixes while skipping optional payloads. Ordinary comments live only in explicitly requested raw strings. Raw Nexus descriptions also need their original translation/command context for equivalent worker parsing; a bare description cannot distinguish TREE from UTREE.

#### Skip burn-in without constructing discarded Trees

```python
import phylotreelib as pt
with pt.Newicktreefile(filecontent="(A[ordinary;comment],B);[&U](C,D);") as tf:
    assert tf.skip_trees(1) == 1
    retained = tf.readtree()
assert retained.leaves == {"C", "D"} and retained.rooted is False
```

`fast_treecount(filename, fileformat="newick")` uses the same lexical framing. Count/skip/raw paths do not decode payloads or validate topology/declaration semantics. Removed `strip_comments` arguments raise TypeError; `readtree(False)` and `get_treestring()` now return intact descriptions. Replace removed `below_root` with structured root-node/root-suffix fields or retained raw text; see the [input API](api.md#input-parsing-and-complete-descriptions).

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
tree = pt.Tree.randtree(ntips=20)
```

#### Random tree with specified leaf names

```python
import phylotreelib as pt
tree = pt.Tree.randtree(leaflist=["Human", "Chimp", "Bonobo", "Gorilla", "Orangutan"])
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

Public `reroot()`, `rootmid()`, `rootout()` and `rootminvar()` establish `tree.rooted=True` on success, including a successful no-op at the existing root. `unroot()` deliberately establishes False. Mechanical `collapse_bifurcating_root()` and deprecated `deroot()` preserve False/None, and clear True to None only when a root is actually removed. Root replacement/removal clears its suffix; a true same-root no-op retains it. Retained edge/node annotations are not automatically recalculated for their new scientific interpretation.

#### Midpoint rooting

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5)
print(tree)
tree.rootmid()
print(tree)
```

#### Minimum-variance rooting

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=5)
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

#### Establish an unrooted interpretation

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
tree.unroot()
print(tree.nexus())  # [&U]; the traversal representation still has a root node.
```

#### Remove one leaf

Pruning preserves the remaining Tree's rootedness, including promotion of the surviving root child. The old suffix is cleared only if that root is removed. Surviving node annotations stay attached; deleted-node data is removed, and serially merged branches keep summed lengths with fresh non-length fields.

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

**Note:** `subtree()` does not modify the source. An extracted Tree inherits its rootedness and independently copies retained mutable annotation data. A proper internal subtree has `root_suffix=None`; its former incoming Branchstruct is returned separately as `basal`. Extraction at the source root independently copies that root suffix and returns `basal=None`. The leaf case remains `(leaf_name, copied_basal_branch)` without a new wrapper. A basal label/stem is never promoted implicitly into an extracted root suffix.

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

This operation inserts a node *on* an existing edge and divides its structural length. The default `copy_attrs="none"` gives both segments fresh non-length fields. Explicit `"lower"` copies annotations to the original-child segment; `"both"` copies them independently to both. Values are retained without scaling their scientific meaning. A supplied `upper_branchstruct` is still installed directly and receives the calculated upper length.

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

`graft()` inserts a point on an existing target edge and attaches a Tree through its root, or a single leaf name. It preserves recipient rootedness and root suffix. Passing a Tree directly consumes its dictionaries, branches and node metadata; use the recipient afterward. Supply an independent copy when preserving the original, as below. Explicit connecting Branchstructs also transfer directly and should be copied first if an independent original is needed.

`connect_branchstruct` takes precedence over `connect_length`; otherwise a fresh connection uses `connect_length`, default zero. Donor suffix labels/stems/annotations do not supply connection data automatically. To reuse a donor suffix deliberately, give it a meaningful numeric length and pass it explicitly as the connecting struct.

```python
import phylotreelib as pt

left = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,E:2);")
right = pt.Tree.from_string("(X:1,(Y:1,Z:1):1);")

# Choose the edge in `left` where the graft point should be inserted:
child = left.find_mrca({"C", "D"})
parent = left.parent(child)

print(left)

# Preserve right by grafting an independent copy; that supplied copy is consumed.
# Place it 30% down the selected edge, with a separate connection of length 0.25.
left.graft(right.copy_treeobject(), parent=parent, child=child,
           frac_from_parent=0.3, connect_length=0.25,
           graftlabel="grafted_")

print(f"right grafted between nodes {parent} and {child}")
print(left)
```

#### SPR with specified prune and regraft points

Complete SPR preserves rootedness. The suffix stays when the original root survives and is cleared if pruning promotes another node. Moved-node data follows final IDs; merged, split and connecting edges have fresh non-length fields, while internal moved edges and unrelated edges retain their annotations. Retention does not prove continued scientific validity. Total length uses the existing basal-length connection, plus the other root-child length when pruning a direct child of a bifurcating root.

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
tree = pt.Tree.randtree(ntips=5)
print(tree.newick(precision=2))
prune_node = next(iter(tree.leaves))
prune_node, regraft_node = tree.spr(prune_node=prune_node)
print(f"prune_node; {prune_node}\tregraft_node = {regraft_node}")
print(tree.newick(precision=2))
```

#### Fully random SPR moves

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=7)
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
tree = pt.Tree.randtree(ntips=200)
tree.prune_maxlen(nkeep=30, keeplist=["s001", "s120", "s178"])
```

#### Ask which leaves would be kept without actually pruning, color those on original tree

```python
import phylotreelib as pt
tree = pt.Tree.randtree(ntips=50)
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
t1 = pt.Tree.randtree(ntips=30)
t2 = t1.copy_treeobject()
t2.spr()
print(t1.treedist_RF(t2))
print(t1.treedist_RF(t2, normalise=True))
```

#### Rooted RF distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=20)
t2 = t1.copy_treeobject()
t2.rootmid()
print(t1.treedist_RF(t2, rooted=False))
print(t1.treedist_RF(t2, rooted=True))
```

#### Path-difference distance

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=30)
t2 = t1.copy_treeobject()
for _ in range(3):
    t2.spr()
print(t1.treedist_pathdiff(t2))
```

#### Request detailed RF output

```python
import phylotreelib as pt
t1 = pt.Tree.randtree(ntips=10)
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

**Note:** For command-line use, the separate tool `sumt` provides the same functionality with a simpler interface. The examples below show the explicit low-level workflow in `phylotreelib` (TreeSummary --> SummaryTreeBuilder --> TreePostProcessor --> Tree.newick() / Tree.nexus()).

* `sumt` GitHub page: [https://github.com/agormp/sumt](https://github.com/agormp/sumt/)

More detailed discussion of summary tree methods is in `docs/api.md`. Input observations use independent accumulators, and constructed summary output Trees have independent branch/annotation data. Retained representatives and consumed partial summaries have the separate contracts below.

#### Merge completed partial summaries

Complete any diagnostics that need separate summaries before merging. A successful `combined.update(part)` combines matching branch, clade, root and topology records and may adopt unseen records directly, including histograms, interval dictionaries and child-clade-pair sets. Use the recipient for subsequent accumulation and statistics; only the donor's scalar count remains safe for progress reporting. No automatic copy, donor reset or reuse guard is added.

```python
import phylotreelib as pt

combined = pt.TreeSummary(trackbips=True, trackblen=True)
part = pt.TreeSummary(trackbips=True, trackblen=True)
combined.add_tree(pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);"))
part.add_tree(pt.Tree.from_string("((A:2,C:3):4,B:5,D:6);"))

combined.update(part)
print("merged trees:", len(part))  # Scalar progress count remains available.
print("combined trees:", len(combined))
summary_branches = combined.bipartsummary
```

Do not add observations to `part`, query its statistics or merge it into another independent result afterward. Recipient updates/finalization may change transferred records.

#### Edit a copy of a retained representative

With clade/bipartition topology tracking and `store_trees=True`, the first Tree for each topology is retained by reference. Leave its topology, names, lengths, rootedness and metadata unchanged while any summary retains it. Read-only queries, printing and normal cache work remain allowed. Copy before editing:

```python
import phylotreelib as pt

sample = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
summary = pt.TreeSummary(trackbips=True, tracktopo=True, store_trees=True)
summary.add_tree(sample)
record = next(iter(summary.biptoposummary.values()))

editable = record.tree.copy_treeobject()
editable.setlength(editable.parent("A"), "A", 2.0)
print(editable.newick())
```

The original sample remains retained. Its no-edit lifetime follows references transferred by summary merges. A Tree may be retained in both topology dictionaries; storage makes no automatic snapshot and adds no mutation guard. The default `store_trees=False` has no representative-storage lifetime restriction. A retained sample's arbitrary annotations are not aggregate summary statistics.

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

# 4) Select labels and metadata for this output
print(contree.newick(
    labelfield="bipartition_cred",
    branch_attributes=("bipartition_cred", "length", "length_sd"),
    print_meta=True, precision=7,
))
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

print(contree.newick(
    labelfield="bipartition_cred",
    branch_attributes=("bipartition_cred", "length", "length_sd"),
    print_meta=True, precision=7,
))
```

#### Build an MCC tree with common-ancestor heights and credible intervals

Both passes must use the same post-burn-in sample of rooted clock trees. The example file below is a demonstration sample; when using your own posterior, discard burn-in consistently in both loops.

```python
# Download an example tree-sample file (run these 3 lines only once)
import urllib.request
URL = "https://raw.githubusercontent.com/agormp/phylotreelib/main/tests/data/primate-mtDNA.trees"
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
plan = pt.CAHeightEstimator.build_plan(mcctree, trackci=True, ci_probs=(0.8, 0.95))
est = pt.CAHeightEstimator(plan, trackci=True)

with pt.Nexustreefile("primate-mtDNA.trees") as tf:
    for tree in tf:
        est.add_tree(tree)

est.write_into(mcctree)
mcctree.set_blens_from_heights()

# Annotate support etc from TreeSummary (does not change heights you just wrote)
tpp = pt.TreePostProcessor(ts)
tpp.annotate_sumtree(mcctree)

# Record this rooted-clock interpretation after construction/processing.
mcctree.rooted = True
# Export all annotations, including root heights, medians, and credible intervals.
with open("mcc_ca_intervals.nexus", "w") as out:
    out.write(mcctree.nexus(print_meta=True))
```

The exported intervals are named `height_80%_CI` and `height_95%_CI`; each holds `{lower,upper}`. They are central, equal-tailed credible intervals, not HPD intervals. `CAHeightEstimator` uses MRCA heights from every sampled tree, including trees where the target clade was not monophyletic. Quantiles are approximate unless a node's height is effectively constant, in which case its mean, median, and interval endpoints coincide.

For a restricted export, select the complete interval key alongside the desired scalar fields:

```python
with open("mcc_ca_heights_only.nexus", "w") as out:
    out.write(mcctree.nexus(
        print_meta=True,
        node_attributes=("height", "height_median", "height_95%_CI"),
        branch_attributes=(),
    ))
```

For clade-height intervals in a single pass, see the [README example](../README.md#compute-credible-intervals-and-export-a-summary-tree). For direct branch-length intervals on a consensus tree, see the [API example](api.md#106-what-the-tracked-quantiles-are-used-for).

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

print(mbctree.newick(
    labelfield="bipartition_cred",
    branch_attributes=("bipartition_cred", "length"),
    print_meta=True, precision=7,
))
```


---

### 10. Exporting with labels, metadata, and color

Pass options directly to each `newick()` or `nexus()` call. Output calls do not store formatting preferences. `tree.nexus(print_meta=True)` includes eligible populated native fields and `extra_meta` entries, including intervals. Supply exact `node_attributes`/`branch_attributes` keys to restrict output, for example `("height", "height_95%_CI")`; use an empty tuple to suppress metadata for that owner. The former `.ci` mapping, `"ci"` group and writer `ci_labels` argument have been removed. Both Tree and TreeSet use full stored precision by default (`precision=None`); pass `precision=7` for seven significant digits. Omit `print_meta` for plain output.

Whole-value None is unset, `pt.MISSING` emits `?`, and valid empty values, zero and False are retained. Use MISSING rather than None inside collections. Branch annotations require printed lengths: with `printdist=False`, pass `branch_attributes=()` or disable metadata. This restriction does not suppress ordinary internal labels. See the [API storage notes](api.md#annotation-storage-and-ordinary-node-records) for native-field mappings, dictionary access and selected-value validation.

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
	tree nexus_tree = ((1:1.0,2:1.0):1.0,(3:1.0,4:1.0):1.0);
end;
```

#### Print branch metadata in BEAST/FigTree-style comments

```python
import phylotreelib as pt
tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
node = tree.find_mrca({"A", "B"})
parent = tree.parent(node)
tree.set_branch_attribute(parent, node, "posterior", 0.98)
print(tree.newick(
    print_meta=True,
    node_attributes=(),
    branch_attributes=("posterior",),
    precision=6,
))
```

Output:

```
((A:1,B:1):1[&posterior=0.98],(C:1,D:1):1);
```

#### Declare rooted or unrooted interpretation

Set the tree's known interpretation before export; its traversal root alone does not establish biological rootedness:

```python
tree.rooted = True
print(tree.nexus())                       # [&R]
tree.rooted = False
print(tree.nexus())                       # [&U]
print(tree.nexus(print_rootedness=False))  # Omit the known declaration
```

`tree.rooted = None` omits the declaration when interpretation is unknown. Nexus TREE leaves rooted/unrooted interpretation unspecified without either declaration. See the [Nexus TREES specification, p. 615](https://www.phylo.org/sub_sections/nexus_document.pdf#page=26) (accessed 2026-10-28). UTREE supplies False. The boolean `print_rootedness` switch defaults to True and works independently of `print_meta`; it never changes state or geometry. `TreeSet.nexus()` uses each tree's own state. The old `nexus(rooted=...)` value argument has been removed.

Constructor and editing rootedness transitions are implemented: clade/MCC/HIPSTR/mrHIPSTR and ancestral-model random trees are rooted, bipartition/consensus/MBC construction is unrooted, and deliberate rooting establishes True. Generic construction and undeclared Newick/TREE input remain None. Input declarations and UTREE context now establish the supported interpretation independently of optional metadata decoding. Root labels/stem lengths and annotations use `root_suffix`; see the [API root-information notes](api.md#root-information-and-copying).

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

### 12. Parsimony

phylotreelib implements Hartigan's (1973) parsimony algorithm, extended to handle
multifurcations and internal nodes with observed states. The implementation follows
[Papamichail et al. (2017)](https://doi.org/10.1016/j.compbiolchem.2017.03.007).

Note: the aim here is to reconstruct ancestral states and count changes on a single given
tree, not to search tree-space for the most parsimonious tree.

Three methods form a pipeline of increasing summarization:

- `parsimony_possible_states()` — optimal state-sets at every node (Hartigan up+down passes)
- `parsimony_assign_fits()` — pick one concrete state per node (one optimal reconstruction)
- `parsimony_count_changes()` — parsimony score + directed change counts

The first two return an annotated copy of the tree; the third returns a score and a
dictionary of change counts. None of them modify the original tree.

States can be single strings (one character in the cladistic sense) or tuples of strings
(multiple independent characters analyzed jointly on the same tree). Single-site and
multi-site examples are shown below.

#### Set up a tree with observed states on leaves

Parsimony methods read `.state` from each node's `Nodestruct` in `tree.nodedict`.
Leaves must have a non-empty `.state`; internal nodes may retain the default None (to be inferred)
or receive an observed state (which the algorithm will respect). An explicitly assigned empty string remains valid annotation data, but does not supply an observed parsimony state.
```python
import phylotreelib as pt

tree = pt.Tree.from_string(
    "((Mink_1:1,Mink_2:1):1,(Human_1:1,Human_2:1):1,(Bat_1:1,Bat_2:1):1);"
)

# Assign observed states to leaves
for leaf in tree.leaves:
    if leaf.startswith("Mink"):
        tree.set_node_attribute(leaf, "state", "mink")
    elif leaf.startswith("Human"):
        tree.set_node_attribute(leaf, "state", "human")
    else:
        tree.set_node_attribute(leaf, "state", "bat")
```

#### Find optimal state-sets at internal nodes
```python
result = tree.parsimony_possible_states()

# result is a copy of tree; inspect the optimal_set at each internal node
for node in result.sorted_intnodes():
    ns = result.nodedict[node]
    print(f"node {node}: optimal_set={ns.optimal_set}  "
          f"primary_set={ns.primary_set}  secondary_set={ns.secondary_set}")
```

`optimal_set` gives the full set of states that a node can have on any maximum-parsimony
reconstruction. `primary_set` and `secondary_set` are the intermediate Hartigan sets.

#### Assign one concrete fit per node
```python
fitted = tree.parsimony_assign_fits()

for node in fitted.nodes:
    ns = fitted.nodedict[node]
    print(f"node {node}: fit={ns.fit}  wasambig={ns.wasambig}")
```

When a node is ambiguous (more than one state is equally parsimonious), the algorithm
picks randomly. Use `fitpref` to prefer a specific state at ambiguous sites:
```python
fitted_prefer_human = tree.parsimony_assign_fits(fitpref="human")

for node in fitted_prefer_human.sorted_intnodes():
    ns = fitted_prefer_human.nodedict[node]
    print(f"node {node}: fit={ns.fit}  wasambig={ns.wasambig}")
```

#### Count changes (parsimony score + direction)
```python
pscore, countdict = tree.parsimony_count_changes()

print(f"Parsimony score: {pscore}")
for (from_state, to_state), count in sorted(countdict.items()):
    print(f"  {from_state} -> {to_state}: {count}")
```

With `fitpref`:
```python
pscore, countdict = tree.parsimony_count_changes(fitpref="human")
print(f"Parsimony score (fitpref='human'): {pscore}")
for (from_state, to_state), count in sorted(countdict.items()):
    print(f"  {from_state} -> {to_state}: {count}")
```

#### Parsimony with an observed state on an internal node

If an internal node has an observed state, set it before calling any parsimony method.
The algorithm will treat it as fixed (like a leaf state).
```python
import phylotreelib as pt

tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")

# Label leaves
for leaf in ["A", "B"]:
    tree.set_node_attribute(leaf, "state", "X")
for leaf in ["C", "D"]:
    tree.set_node_attribute(leaf, "state", "Y")

# Also fix the root's state
tree.set_node_attribute(tree.root, "state", "X")

pscore, countdict = tree.parsimony_count_changes()
print(f"Parsimony score: {pscore}")
for (f, t), c in sorted(countdict.items()):
    print(f"  {f} -> {t}: {c}")
```

#### Multi-site parsimony

When each leaf carries observations at multiple independent sites (characters), you can
pass the state as a tuple of strings instead of a single string. The Hartigan algorithm
then runs independently at each site in a single tree traversal, and the parsimony score
and change counts are aggregated across all sites.

All labelled nodes must use tuples of the same length. Unlabelled nodes use None by default.
An explicit empty string `""` is also treated as unlabelled.

```python
import phylotreelib as pt

tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")

# Two sites: a "host" character and a nucleotide character
tree.set_node_attribute("A", "state", ("mink", "G"))
tree.set_node_attribute("B", "state", ("mink", "A"))
tree.set_node_attribute("C", "state", ("human", "G"))
tree.set_node_attribute("D", "state", ("human", "A"))

pscore, countdict = tree.parsimony_count_changes()
print(f"Parsimony score (summed over 2 sites): {pscore}")
for (f, t), c in sorted(countdict.items()):
    print(f"  {f} -> {t}: {c}")
```

Note: `countdict` keys are individual state strings (e.g., `("mink", "human")`), not
tuples of tuples. Changes at different sites that happen to involve the same pair of
states are summed together.

#### Inspect per-node results in multi-site mode

When states are tuples, the node attributes on the returned tree are also tuple- or
list-shaped:
```python
import phylotreelib as pt

tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
tree.set_node_attribute("A", "state", ("X", "1"))
tree.set_node_attribute("B", "state", ("X", "2"))
tree.set_node_attribute("C", "state", ("Y", "1"))
tree.set_node_attribute("D", "state", ("Y", "2"))

result = tree.parsimony_possible_states()

for node in result.sorted_intnodes():
    ns = result.nodedict[node]
    # optimal_set is a list of sets, one per site
    print(f"node {node}: optimal_set={ns.optimal_set}")

fitted = tree.parsimony_assign_fits()

for node in fitted.sorted_intnodes():
    ns = fitted.nodedict[node]
    # fit is a tuple of strings, wasambig is a tuple of bools
    print(f"node {node}: fit={ns.fit}  wasambig={ns.wasambig}")
```

#### Multi-site fitpref

`fitpref` accepts either a single string (broadcast to all sites) or a tuple with one
entry per site. Per-site entries may be `None` for no preference at that site.
```python
import phylotreelib as pt

tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);")
tree.set_node_attribute("A", "state", ("X", "1"))
tree.set_node_attribute("B", "state", ("Y", "1"))
tree.set_node_attribute("C", "state", ("X", "2"))
tree.set_node_attribute("D", "state", ("Y", "2"))
tree.set_node_attribute("E", "state", ("X", "1"))
tree.set_node_attribute("F", "state", ("Y", "2"))

# Broadcast: prefer "X" at every site where it is in the optimal set
pscore, countdict = tree.parsimony_count_changes(fitpref="X")
print(f"Parsimony score (fitpref='X'): {pscore}")

# Per-site: prefer "X" at site 0, no preference at site 1
pscore, countdict = tree.parsimony_count_changes(fitpref=("X", None))
print(f"Parsimony score (fitpref=('X', None)): {pscore}")
```

Note: `fitpref` only affects the choice among equally parsimonious reconstructions. It
does not change the parsimony score when the reconstruction is unambiguous.



### Editing lengths and retained annotations

`setlength()` and `set_branch_attribute(..., "length", value)` change structural lengths and refresh derived length caches without clearing stored node/branch annotations. Caller scaling loops follow the same rule; `set_blens_from_heights()` uses explicitly stored heights and rejects missing required values. An imported annotation named `length` stays distinct from a colon length, and stored node heights stay distinct from `nodeheightdict` derived from branch lengths.

For low-level batch writes directly to `tree.child_dict[parent][child].length`, call `tree.clear_length_caches()` after the batch. Ordinary metadata assignments do not invalidate topology caches or silently recalibrate statistics. Reassess retained scientific annotations explicitly after topology or length edits.

`Branchstruct.merge(other)` creates a fresh length-only result, and its unused `check_compat` argument has been removed. `remove_branch()` redistributes the removed structural length equally onto fresh descendant-edge records, retaining unrelated edges and surviving-node data. Existing edges that only reverse during rerooting keep their labels/annotations, whose directional or rooted-clade meaning may need reassessment.

### 13. Upgrading existing 1.x scripts

These examples describe the changes introduced in 2.0. For the 3.0 changes, including removal of stored printing settings, read the [3.0 migration guide](api.md#upgrading-to-30).

phylotreelib 2.0 introduced a few intentional API breaks to make tree editing, summary trees, and printing more explicit and less "magic".

If you maintain scripts that call phylotreelib directly, the sections below show the key signature/return changes and how to update call sites.

#### Terminology rename: `depth` → `height`

In v1.x, the distance from the tips to an internal node was called **depth**. In v2.0.0 this has been renamed to **height** throughout the library. The rename is purely mechanical: a search-and-replace across your scripts using the table below is sufficient.

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

#### Summary trees

`TreeSummary.compute_sumtree(...)` is still available as a convenience wrapper (so many 1.x workflows continue to work), but the preferred style is now the explicit pipeline:

- `TreeSummary` --> `SummaryTreeBuilder` --> `TreePostProcessor` (+ optional `CAHeightEstimator`) --> `Tree.newick()` / `Tree.nexus()`

See the [README summary example](../README.md#summarize-posterior-tree-samples) and the summary-tree recipes above.

#### Tree.graft: signature change

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

#### Tree.spr: signature change (and length preservation)

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

#### Tree.subtree / Tree.prune_subtree: return value change

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

#### Tree.reroot: signature change

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

#### Root-collapsing helper rename

**Old (1.x):**

```python
Tree.deroot(self)
```

**New (2.x):**

```python
Tree.collapse_bifurcating_root(self)
```

(And `Tree.deroot()` is kept as a deprecated alias.)

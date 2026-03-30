## phylotreelib API reference and advanced notes

This file is a curated API reference for the public interface to `phylotreelib`, followed by a few advanced notes on summary-tree computation and quantile tracking.

It is intentionally selective: the goal is to document the classes and functions that typical users should rely on directly.

---

### 1. Public API overview

The main public entry points are:

- `Tree`
- `Treefile`, `Newicktreefile`, `Nexustreefile`
- `TreeSet`
- `TreeSummary`
- `Distmatrix`
- `TreeError`
- `build_sumtree`
- `configure_basic_printing`
- `configure_sumtree_printing`

Helper classes such as `Bipartition`, `Clade`, `Branchstruct`, `Nodestruct`, `PrintSpec`, `SummaryTreeBuilder`, `TreePostProcessor`, `CAHeightEstimator`, and `QuantileAccumulator` are also part of the module, but most users do not need them for ordinary scripting.

---

### 2. Tree input and output

#### `Treefile(filename, strip_comments=True)`

Factory that autodetects file format and returns either a `Newicktreefile` or `Nexustreefile`.

Use this when you want one entry point regardless of the underlying tree-file format.

#### `Newicktreefile(filename=None, filecontent=None, strip_comments=True)`

Class representing a Newick tree file.

#### Building summary trees (preferred)

The summary-tree machinery is split into three roles:

- **TreeSummary**: accumulates frequencies and (optionally) branch length / node height statistics from many trees.
- **SummaryTreeBuilder**: chooses a *summary topology* (consensus, MCC, MBC, HIPSTR).
- **TreePostProcessor**: annotates a chosen tree with support values, length/height statistics, and (optionally) root credibility.

For most users, the recommended entry point is the top-level helper:

#### `phylotreelib.build_sumtree(...)`

```python
import phylotreelib as pt

# 1) Accumulate summaries
ts = pt.TreeSummary(trackbips=True, trackblen=True)
with pt.Nexustreefile("posterior.trees") as tf:
    for tree in tf:
        ts.add_tree(tree)

# 2) Build + annotate the summary tree (topology + rooting + lengths/heights + support)
sumtree = pt.build_sumtree(ts, treetype="con", blen="biplen", rooting="mid")

# 3) Configure output (optional)
pt.configure_sumtree_printing(sumtree, treetype="con", blen="biplen", print_meta=True)

print(sumtree.newick())
```

Parameters
- `treetype`: `"con"`, `"all"`, `"mcc"`, `"mbc"`, `"hip"`, `"mrhip"`
- `blen`: `"biplen"`, `"cladeheight"`, `"caheight"`, `"input"`, `"none"`
- `rooting`: `"mid"`, `"minvar"`, `"og"`, or `None` (meaning: keep as-is / leave unrooted)

Notes
- `"caheight"` needs access to the original tree sample again. Pass `count_burnin_filename_list`
  (the same structure used by `sumt`: tuples of `(count, burnin, filename)`), or do a second pass
  yourself with `CAHeightEstimator`.
- Printing is intentionally **not** configured by `build_sumtree`; use `Tree.set_print_spec()`,
  `configure_sumtree_printing()`, or `configure_basic_printing()`.

**v2.0.0 note:** The `blen` values `"meandepth"` and `"cadepth"` from v1.x have been renamed to `"cladeheight"` and `"caheight"` respectively.

#### Lower-level building blocks

If you want more control (or are writing another CLI like `sumt`), you can call the building blocks directly:

```python
stb = pt.SummaryTreeBuilder(ts)
tpp = pt.TreePostProcessor(ts)

tree = stb.contree(cutoff=0.5, allcompat=False)   # topology
tree.rootmid()                                    # rooting (Tree method)
tree = tpp.annotate_sumtree(tree)                 # support + length/height stats
```

#### Backwards-compatible wrappers on `TreeSummary`

`TreeSummary` still exposes convenience wrappers like `ts.contree()`, `ts.max_clade_cred_tree()`,
and `ts.annotate_sumtree()`. These exist for compatibility with older scripts, but newer code
should prefer `SummaryTreeBuilder` / `TreePostProcessor` (or simply `build_sumtree`).

#### Summary properties

- `bipartsummary`
- `cladesummary`
- `rootbipsummary`
- `biptoposummary`
- `cladetoposummary`
- `sorted_biplist`
- `sorted_rootbips`

---

### 6. `Distmatrix`

Distance matrix class for distance-based tree reconstruction.

#### Constructors

- `Distmatrix.from_distdict(distdict)`
- `Distmatrix.from_distfile(distfilename)`
- `Distmatrix.from_numpy_array(nparray, namelist)`

#### Main methods

- `nj()` — neighbor joining
- `upgma()` — UPGMA
- `getdist(name1, name2)`
- `setdist(name1, name2, dist)`
- `rename(oldname, newname)`
- `clean_names(illegal=',:;()[]', rep='_')`
- `avdist()`

---

### 7. Top-level helpers

#### `build_sumtree(treesummary, treetype='con', blen='biplen', rooting=None, og=None, count_burnin_filename_list=None)`

Canonical top-level function for building and annotating a summary tree from a populated `TreeSummary`.

Conceptually, this function:

1. chooses a summary-tree topology
2. roots it if requested
3. assigns branch lengths and/or node heights
4. annotates support and associated statistics

#### `configure_basic_printing(...)`

Convenience function for setting a tree's `PrintSpec`.

Use this to choose:

- floating-point precision
- whether to print lengths and labels
- which branch attribute acts as the label
- whether to emit metadata comments
- which node and branch attributes to include in the comments

#### `configure_sumtree_printing(...)`

Convenience function that chooses sensible printing defaults for summary trees, depending on summary-tree type and branch-length convention.

---

### 8. Exceptions

#### `TreeError`

Exception raised for tree-related problems.

Typical uses include parse errors, invalid node references, incompatible operations, or illegal constructor arguments.

---

### 9. Advanced notes: summary-tree methods

This section gives a high-level conceptual explanation of the main summary-tree methods.

#### 9.1 Consensus trees

A consensus tree is built from splits observed across many trees.

In the majority-rule version, only bipartitions with posterior frequency above the chosen cutoff are retained. This gives a direct summary of commonly supported branches.

If `allcompat=True`, additional compatible bipartitions can be added even if they fall below the cutoff, producing a more resolved tree without violating compatibility.

Use consensus trees when your main goal is to summarize branch support rather than to select one specific sampled topology.

#### 9.2 MCC and MBC trees

##### Maximum clade credibility (MCC)

The MCC tree is the sampled tree whose product of clade posterior frequencies is maximal, or equivalently whose sum of log clade frequencies is maximal.

This approach prioritizes trees built from highly supported clades.

##### Maximum bipartition credibility (MBC)

The MBC tree is analogous but uses branch bipartitions instead of rooted clades. It therefore evaluates trees in terms of branch support rather than rooted clade support.

##### Practical difference

- MCC is usually the natural choice for rooted posterior trees, especially clock trees.
- MBC can be useful when input trees are unrooted or bipartition support is the quantity of interest.

Both methods choose one sampled topology rather than constructing an abstract consensus topology.

#### 9.3 HIPSTR and mrHIPSTR

`phylotreelib` also supports HIPSTR-style summary trees, following the idea of highest independent posterior subtree reconstruction.

The key idea is to combine local posterior subtree information in a principled way to construct a summary topology.

Practical rule of thumb:

- use HIPSTR when you want a summary-tree method based on posterior subtree reconstruction rather than a single sampled tree's global score
- use mrHIPSTR for the majority-rule variant

#### 9.4 Root summaries

When `trackroot=True`, `TreeSummary` can record where the root occurred across sampled trees.

This allows:

- placing the root of a summary tree on the most frequently observed root bipartition
- annotating each branch with `rootcred`, the posterior frequency with which that branch carried the root

This is useful when the posterior distribution contains uncertainty about root placement.

#### 9.5 Branch lengths versus node heights

Summary trees can be annotated using different conventions.

##### Bipartition-based branch lengths

A branch in the summary tree corresponds to a bipartition, and branch length is summarized directly across all sampled trees where that bipartition occurred.

This is natural for consensus trees. Use `blen="biplen"` and `trackblen=True` in `TreeSummary`.

##### Node-height-based branch lengths

For rooted clock trees, one often wants node heights rather than direct bipartition branch lengths. Node height is the distance from the tips to the node — i.e., "time before most recent leaf".

Two closely related conventions appear in the package:

- **Mean clade heights** (`blen="cladeheight"`, requires `trackheight=True`): branch lengths are derived from the mean height of each clade, computed across all sampled trees where that clade was monophyletic.
- **Common-ancestor heights** (`blen="caheight"`): the height of a clade on the summary tree is estimated from the MRCA of those leaves in *every* sampled tree, even when the clade is not monophyletic in a particular sample. This corresponds to the common-ancestor-height convention used by BEAST's tree annotation workflows.

**v2.0.0 note:** These options were called `blen="meandepth"` and `blen="cadepth"` in v1.x. The parameter in `TreeSummary` that enables node-height tracking was called `trackdepth`; it is now `trackheight`.

#### 9.6  Branch labels are treated as branch attributes

Tree objects created from string representations (newick or nexus) will interpret labels as belonging to the branch (bipartition) not a node at either end of that branch. This affiliation will remain after rerooting.


---

## 10. Advanced notes: quantile tracking and the QuantileAccumulator

When `trackci=True`, `TreeSummary` tracks approximate quantiles for branch lengths and/or node heights.

This is done using the `QuantileAccumulator` class.

### 10.1 Why an approximate quantile accumulator?

For large posterior tree samples, storing every observed branch length or node height for every summary feature can become expensive in both memory and merge cost.

The package instead uses a mergeable approximate summary so that it can:

- process data in one pass
- keep memory use bounded
- merge partial summaries efficiently, which is helpful for large workflows and possible multiprocessing

### 10.2 Basic idea

The accumulator uses a log-bucket histogram.

For each positive observed value `x`:

1. write `x = m * 2^e`, where `m` is the mantissa and `e` the exponent
2. discretize the mantissa into one of a fixed number of sub-buckets
3. use the pair `(e, mantissa bucket)` as the bucket key
4. store only counts per bucket

Because the bucket keys are ordered in the same way as the underlying values, one can recover approximate quantiles by walking through buckets in order and locating the bucket that contains the desired rank.

### 10.3 Why logarithmic bucketing?

Logarithmic bucketing gives approximately constant **relative** resolution across scales.

That is often a good fit for branch lengths and node heights, which may span multiple orders of magnitude.

### 10.4 Mergeability

Two accumulators can be merged simply by adding bucket counts. This is important because it means approximate quantiles remain easy to combine across partial summaries.

### 10.5 Accuracy trade-off

The method is approximate, not exact.

The benefit is that it avoids storing all raw values while still providing useful credible intervals and quantiles for summary-tree annotation.

If exact quantiles are needed for a particular analysis, they should be computed from the raw values outside the accumulator framework.

### 10.6 What the tracked quantiles are used for

Quantile summaries are typically used to populate metadata such as:

- branch-length credible intervals
- node-height credible intervals
- medians and related robust summaries

These can then be printed into Newick/NEXUS metadata comments via `configure_basic_printing()` or `configure_sumtree_printing()`.

---

## 11. Public vs internal API

For everyday use, it is reasonable to treat the following as the stable public interface:

- `Tree`
- `Treefile`, `Newicktreefile`, `Nexustreefile`
- `TreeSet`
- `TreeSummary`
- `Distmatrix`
- `TreeError`
- `build_sumtree`
- `configure_basic_printing`
- `configure_sumtree_printing`

The remaining classes are useful, but many of them are better regarded as advanced helpers or implementation-level infrastructure unless they are explicitly needed by your workflow.

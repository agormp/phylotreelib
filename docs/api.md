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
- `MISSING`
- `build_sumtree`

Helper classes such as `Bipartition`, `Clade`, `Branchstruct`, `Nodestruct`, `SummaryTreeBuilder`, `TreePostProcessor`, `CAHeightEstimator`, and `QuantileAccumulator` are also part of the module, but most users do not need them for ordinary scripting.

#### Upgrading to 3.0

Version 3.0.0 changes the public API and requires updating some 2.x callers. The linked sections describe the full behavior.

| Earlier use | Replacement or changed behavior |
|---|---|
| `PrintSpec`, `configure_*_printing()`, stored print-setting methods | Pass options directly to each `newick()` or `nexus()` call. |
| Public `.ci`, writer `ci_labels`, `"ci"` output group | Store intervals in `extra_meta` under complete keys such as `height_95%_CI`, and select exact keys. Literal `ci` remains ordinary metadata. |
| Reader/factory `strip_comments` | Use keyword-only `parse_meta=True` for decoded annotations, or raw reading for original text. |
| Reader/Tree `below_root` | Use root-node data in `nodedict[root]`, root-label/stem data in `root_suffix`, or retain raw text explicitly. |
| `nexus(rooted=...)` | Set `tree.rooted`, then use `print_rootedness` to include or suppress its declaration. False suppresses output, rather than declaring unrooted. |
| Omitted printing precision | Full stored precision is now the default. Pass `precision=7` for seven significant digits. |
| Positional metadata/declaration switches | Pass `print_meta` and `print_rootedness` by keyword. Removed interval-filter calls raise TypeError. |
| `Nodestruct.nleaves` | Use summary records or topology-derived descendant counts. |
| Branch subdivision copying attributes by default | The default is `copy_attrs="none"`. Request `"lower"` or `"both"` explicitly when needed. |
| `Branchstruct.merge(check_compat=...)` | Remove the unused argument. Merging retains summed length and fresh non-length fields. |
| `copy_treeobject(copyattr=False)` retaining node data | Node metadata is now omitted too. |

Before, in 2.x:

```python
pt.configure_basic_printing(tree, precision=7, print_meta=True, branch_attrs=("posterior",))
print(tree.nexus())
```

Now:

```python
print(tree.nexus(precision=7, print_meta=True,
                 node_attributes=(), branch_attributes=("posterior",)))
```

The sections below cover [input/output](#2-tree-input-and-output), [editing](#3-tree-construction-and-editing) and [summary ownership](#5-summary-ownership). Graft donors and merged partial summaries transfer their records. Retained input representatives must remain unchanged while a summary uses them.

For older changes, see the [1.x upgrade recipes](recipes.md#13-upgrading-existing-1x-scripts).

---

### 2. Tree input and output

#### `Tree.newick(...)` and `Tree.nexus(...)`

Printing is controlled by each call. No formatting settings are stored on a tree.

```python
Tree.newick(self, printdist=True, printlabels=True, labelfield="label", precision=None,
            transdict=None, node_attributes=None, branch_attributes=None, *,
            print_meta=False)
Tree.nexus(self, printdist=True, printlabels=True, labelfield="label", precision=None,
           translateblock=False, node_attributes=None, branch_attributes=None,
           colorlist=None, colorfg="#0000FF", colorbg="#000000", *,
           print_meta=False, print_rootedness=True)
```

- `printdist` and `printlabels` control structural branch lengths and internal labels; `labelfield` names the branch field used for labels.
- Omitted precision or `precision=None` uses the shortest decimal representation recovering each supported finite stored float. Explicit numeric precision, such as `precision=7`, retains significant-digit formatting. Annotation floats retain a decimal point or exponent when formatted integrally; integers remain exact. This preserves stored precision and types without recovering input digits or increasing statistical precision.
- `print_meta=True` enables optional metadata. Node and branch selections are independent: None includes eligible populated native fields and `extra_meta` entries, an empty sequence selects none, and an explicit sequence selects exact annotation keys. For example, `node_attributes=("height", "height_95%_CI")` selects those two annotations. Automatic native-field output excludes branch accumulators (`mean`, `M2`, `n`, `quantiles`), temporary heights (`parent_height`, `kid_height`) and the `extra_meta` container; matching dictionary entries remain eligible.
- Whole-value None is unset. Selected `MISSING` emits unquoted `?`; valid empty strings/collections, zero and False remain present. None inside a collection is invalid; use MISSING for unavailable elements. Strings are double-quoted with embedded double quotes doubled, booleans use `true`/`false`, and ordered collections retain their order within braces. Sets export in deterministic order. Square brackets must balance independently within each selected string or quoted key.
- Emittable branch annotations require printed branch lengths. With `printdist=False`, select `branch_attributes=()` or disable metadata. Node annotations and ordinary internal labels remain independently controlled.
- Existing scalar printing settings accept explicit None as their fixed default; `print_rootedness` requires a boolean. Nexus output wraps Newick and adds translate-block, taxon-color and rootedness controls.

`TreeSet.newick()` and `TreeSet.nexus()` forward the same controls and full-precision default to each tree. Nexus uses each tree's own rootedness state.

**3.0 migration:** `PrintSpec`, the `configure_*_printing()` helpers and Tree's stored-print-setting methods have been removed. Pass options directly; former helper `node_attrs`/`branch_attrs` become `node_attributes`/`branch_attributes`. Public `.ci`, the `"ci"` interval-group shortcut and writer `ci_labels` are also removed. Select complete interval keys instead; a literal `ci` annotation selects only itself. Parameters preceding the removed interval-filter slot retain their positions; the following `print_meta`/`print_rootedness` switches are keyword-only. Old keyword or positional interval filters raise TypeError. See the [migration table](#upgrading-to-30).

#### Annotation storage and ordinary node records

Recognized node scalar keys (`height`, `height_sd`, `height_median`, `clade_cred`) use their native fields. Recognized incoming-branch scalar keys (`length_median`, `length_sd`, `posterior`, `freq`, `bipartition_cred`, `rootcred`) do likewise. These keys cannot have competing copies in the same owner's `extra_meta`.

All other annotations, including intervals, live in optional `extra_meta`, initially None, on `Nodestruct`, `Branchstruct` and public `SummaryNodestruct` result records. Complete keys are preserved without inferred prefixes. For example:

```python
tree.set_node_attribute("A", "extra_meta", {
    "height_95%_CI": (1.0, 2.0), "location": "Denmark",
})
tree.nodedict["A"].extra_meta["location"] = "Norway"
```

Native getters/setters read or replace whole fields; they do not redirect arbitrary attribute names into dictionary entries or implicitly merge maps. For output, native scalar fields are authoritative; otherwise a non-None extra entry takes precedence over eligible native-field fallback. A literal entry named `extra_meta` is valid and does not serialize the container itself. Only selected values are validated.

Recognized interval keys match `height_[0-9]+(?:\.[0-9]+)?%_(?:CI|HPD)` on nodes or `length_[0-9]+(?:\.[0-9]+)?%_(?:CI|HPD)` on branches. They accept whole MISSING or two ordered numeric/MISSING endpoints. Numeric scalars and endpoints exclude booleans. Wrong-owner and other keys remain generic; pairs alone do not establish interval semantics. Validation does not impose probability ranges, endpoint ordering or scientific consistency.

`Nodestruct(height=None)` no longer accepts or exposes `nleaves`. Descendant counts remain in `SummaryNodestruct` and topology-derived `Tree.iter_cladeinfo()` results. An explicit annotation named `nleaves` remains generic data; an annotation named `length` does not replace a structural colon length. Stored heights likewise remain separate from branch-derived height caches.

The shared public `MISSING` marker retains identity through copying and pickle. `set_blens_from_heights()` rejects absent or explicitly missing required heights. Fresh `Nodestruct.state` and unestimated `Branchstruct.length_sd` are None; an assigned empty string or zero is distinct data.

#### Root information and copying

`tree.root` is the traversal-node identifier. Set `tree.rooted` to True, False or None to record the interpretation independently of that identifier. `tree.nexus()` defaults to `print_rootedness=True`, emitting `[&R]`/`[&U]` for known state and omitting a declaration for None. `print_rootedness=False` omits a known declaration; it does not declare the tree unrooted. TreeSet uses each tree's own state; neither Newick method has a declaration switch. Printing changes neither interpretation nor geometry, and `print_meta=False` does not suppress a declaration. This replaces the old `nexus(rooted=...)` value argument.

`tree.root_suffix` is None or a Branchstruct for root-label, optional stem-length and stem-annotation data outside the represented topology. Programmatically populated suffix data use `printlabels`/`labelfield`, `printdist` and exact `branch_attributes` selections. Root-node annotations remain in `nodedict[root]`, before a root colon. A label-only suffix must have `length=None`; Branchstruct's ordinary default is zero. Selected stem annotations require an actually emitted stem length, including when a native label would otherwise be included as metadata. Selecting no branch metadata still permits an ordinary root label without a stem.

`copy_treeobject(copylengths=True, copyattr=True)` preserves rootedness under either flag and independently copies retained mutable annotation containers. `copyattr=False` omits node metadata and non-length branch/suffix fields. `copylengths=False` resets a present stem to zero while preserving an absent stem as None; a suffix with no retained information is omitted. This does not establish independent copying of unrelated accumulators.

#### Input parsing and complete descriptions

```python
Tree.from_string(orig_treestring, transdict=None, format="newick",
                 label_attr_name="label", label_type=str, *, parse_meta=False)
Treefile(filename, *, parse_meta=False)
Newicktreefile(filename=None, filecontent=None, *, parse_meta=False)
Nexustreefile(filename=None, filecontent=None, *, parse_meta=False)
```

`Treefile` detects the file format. Concrete readers require exactly one of `filename` or `filecontent`. Normal iteration and `readtree()` return Trees; `readtree()` returns None at clean EOF. Each reader reuses its parser.

`parse_meta=True` decodes key/value annotations at their positional owner. Before a colon they belong to the preceding leaf/internal/root node; after it they belong to that node's incoming branch, or the root's stem. Root labels use the existing configured label field/conversion. Root-node data alone creates no `root_suffix`; a root label without a supplied stem has `length=None`. Structural suffix input is retained in both metadata modes.

Values support exact integers, decimal/scientific floats, lowercase booleans, strings, `?` as MISSING and recursively nested brace collections decoded as tuples. Single/double quoted keys and strings use doubled matching quotes; backslashes remain literal. Quoted scalar values remain strings. Bare tokens are classified completely; `nan`, `inf`, `None` and malformed numeric-looking tokens remain generic strings and fail where a recognized numeric field is required. Square brackets must balance independently within every quoted string/key. Native mappings and interval validation follow the storage rules above.

Disjoint blocks on the same owner merge. Repeated exact decoded keys on that owner raise TreeError, including identical values; node and branch owners remain separate. With metadata enabled, malformed payloads, unsupported directives such as `[&W 0.5]`/NHX and leading key/value annotations without a supported owner raise TreeError. The parser does not infer weights, decode producer-specific dialects or broaden quoted taxon/tree-name/translation syntax.

Both metadata modes interpret one standalone leading `[&R]` or `[&U]`, with case/whitespace tolerance. Repeated or misplaced declarations raise TreeError. `[&R=1]` is an ordinary assignment, never a declaration. A Nexus TREE command without a flag leaves None; UTREE establishes False and rejects a contradictory R flag. Interpretation comes from syntax/context rather than tree shape. Nexus itself leaves an undeclared TREE's interpretation unspecified. See the [Nexus TREES specification, p. 615](https://www.phylo.org/sub_sections/nexus_document.pdf#page=26) (accessed 2026-10-28).

With metadata disabled, optional payloads are skipped without decoding or annotation maps. Ordinary comments are discarded from parsed Trees in both modes. Nexus names/pre-`=` comments and unrelated application blocks are envelope material and are not stored. For example, a header `[p=0.2]` neither supplies node metadata nor weights the observation.

`read_tree_text()` returns the complete original description, including its whitespace, ordinary/ampersand comments, declarations and root suffix, or None at clean EOF. Nexus results start immediately after `=` and exclude the header/block envelope. `readtree(returntree=False)` and `get_treestring()` are thin wrappers with the same result. Raw retrieval frames boundaries without validating topology, values or declarations; raw strings alone do not carry UTREE or translation context. Internal worker transport captures stable command/parser context for each corresponding description.

`skip_trees(n)` skips up to a nonnegative number of descriptions and returns the actual count. `fast_treecount(filename, fileformat="nexus")` counts descriptions using the same lexical framing. These paths construct neither complete descriptions nor Trees, and do not validate annotation/topology semantics. All paths reject unfinished lexical input and unterminated/stray comment brackets; brackets nest even inside quotes in comments.

**3.0 input migration:** removed reader/factory `strip_comments` keywords and positional slots raise TypeError. Use `parse_meta=True` for decoded annotations and raw reading for original description text. Default reads now retain supported declarations. Raw wrappers retain their call/EOF shapes but deliberately preserve comments and suffixes. Removed reader/Tree `below_root` has no alias: use root-node data in `nodedict[root]`, root labels/stems in `root_suffix`, or explicitly retain raw text outside the Tree. Standalone `remove_comments()` retains its existing generic behavior and discards declarations when called explicitly.

For an unchanged supported Tree, selected output and metadata-enabled reimport preserve topology/root placement, selected structural information, emitted annotation owners/keys/types/values, and finite floats at omitted/None precision. Match translation and label settings and relevant command context. Canonical Nexus declarations preserve rootedness; plain Newick or suppressed declarations intentionally omit it. This does not promise the original imported key set, Python container identity, comments/names/header history or source formatting. Native fallback may add selected fields that were never imported; scientific validity after edits is the caller's responsibility.


---

### 3. Tree construction and editing

`Tree.intnodes` returns a live view of internal node IDs in dictionary order. The view is not stored in the Tree.

#### Constructor interpretation

Generic `Tree()`, `from_leaves()` and `from_branchinfo()` leave rootedness None. `from_biplist()`, `from_topology()` and neighbor joining establish False. `from_cladedict()`, UPGMA and ancestral random Yule/coalescent generation establish True, including rate-varying random lengths. Consensus/all-compatible consensus/MBC inherit the unrooted construction; MCC/HIPSTR/mrHIPSTR inherit the oriented clade construction. Deliberate subsequent rooting establishes True. Estimating heights, lengths or support preserves interpretation rather than inferring one, and retained representative Trees keep their original state. Undeclared Newick/TREE input remains unspecified; UTREE or a supported declaration supplies interpretation.

#### Rooting, unrooting and structural root collapse

A successful public `reroot()` selects a rooted interpretation, including when selecting the current root. It retains the existing geometry and does not automatically collapse the previous bifurcating root. Actual root replacement clears the old `root_suffix`; a genuine no-op retains it. Surviving nodes retain their node data, and a new root node starts with normal empty data. Internal computational movement uses the same mechanics without selecting rootedness: displacement clears True to None and preserves False/None; a no-op preserves all states. `match_nodes()` uses this mechanical route on its temporary copies.

`unroot()` establishes False after successful structural preparation. A multifurcating or two-leaf traversal root remains in place, preserving its suffix and node data. `collapse_bifurcating_root()` is mechanical: actual removal clears the suffix and deleted-root node data and changes True to None; False/None remain unchanged. A non-bifurcating no-op preserves all states. Direct structural collapse of a two-leaf tree still raises TreeError. The deprecated `deroot()` alias continues to mean structural collapse; use `unroot()` to declare an unrooted interpretation.

Root selectors (`rootmid()`, `rootminvar()`, `rootout()` and successful `TreePostProcessor.root_maxfreq()`) establish True, including success paths that keep the current root. A genuine no-op retains its suffix; an actual remove/recreate step clears it even if the final geometric position is equivalent. `root_maxfreq()` retains its existing root-position calculation; its unfinished length-fraction/rootcred work remains separate.

#### Branch merging

`Branchstruct.merge(other)` returns a fresh branch whose structural length is the sum of both lengths. Labels, statistics, intervals, extras and accumulator state start from normal defaults; neither input changes. The unused `check_compat` argument has been removed. Root-edge combination and the serial join after removing a non-root bifurcating parent follow this same length-only rule. Node data stays with surviving nodes; deleted-node data is removed. Callers may explicitly assign justified annotations to the result.

#### Branch subdivision

`add_node_on_branch(..., copy_attrs="none")` divides structural length and starts both segments with fresh non-length fields by default. `copy_attrs="lower"` explicitly copies annotations to the original-child segment; `"both"` copies them independently to both segments. Copied statistics are not rescaled. `upper_branchstruct` still installs the supplied object directly and sets its length to the calculated upper length. Automatic graft and reroot subdivisions also start without source annotations. Endpoint-node data remains attached to its nodes, and the inserted node starts with normal empty data. `split_off_children()` moves existing edges and is outside this subdivision rule.

#### Branch removal and length redistribution

`remove_branch()` removes an internal edge and its child node, adding an equal share of the removed structural length to each reattached descendant edge. Total tree length is conserved; pairwise distances generally change. The reattached edges have new endpoints and fresh non-length fields even when the removed length is zero. The deleted node's metadata is removed. Unrelated edges and surviving-node data are retained; no statistic or annotation named `length` is numerically redistributed.

#### Existing-edge reversal

Reversing an existing edge retains its Branchstruct and annotation values unchanged. Its endpoints, length and unrooted split remain, so reversal does not copy or reinterpret annotations. Directional values and rooted-clade support labels may no longer describe the descending side; callers must reassess/recompute them where needed. This attachment rule does not restore data already discarded by a separate merge or subdivision, nor does it guarantee scientific validity after rerooting.

#### Pruning and extraction

Iterating over a Tree yields proper internal subtrees in `sorted_intnodes()` order. Each result is a new Tree with its copied incoming branch in `.basalbranch`. The source root and individual leaves are excluded, so a star tree yields no subtrees. Use `subtree(leaf)` to extract an individual leaf and its basal branch.

Successful `remove_leaf()`, `remove_leaves()`, `prune_subtree()` and pruning wrappers preserve the remaining Tree's rootedness, including surviving-child promotion. This retains the oriented surviving clade and differs from structural root collapse retaining both sides. Actual removal of the old root clears its entire suffix and deleted-node data; surviving-node annotations stay attached. A selection-only wrapper does not mutate state or suffix. Extracted Tree state is inherited before source mutation. Existing valid-retained-Tree geometry and singleton/empty boundaries are unchanged.

A Tree returned by `subtree()` inherits True, False or None from the source at extraction time, including extraction at the source root. `prune_subtree()` obtains that result before editing the source. Selecting a traversal subtree does not infer biological rootedness. Leaf extraction still returns a name string plus its copied basal branch; it creates no Tree state or new wrapper.

`subtree(source.root)` independently copies retained source suffix data and returns `basal_branch=None`. A proper internal subtree has `root_suffix=None` and returns its former incoming branch separately, including its label and annotations. Leaf extraction keeps the same string/copied-basal-branch interface. A basal branch never becomes an implicit root label or stem. Extracted node data, internal edges, the basal branch and a copied suffix isolate supported mutable annotations. `prune_subtree()` obtains the extracted result before source mutation. Pruning the entire source root is outside the existing valid-retained-Tree contract.

#### Grafting and donor ownership

A successful `recipient.graft(donor, ...)` consumes donor contents: ID/name renaming and direct adoption of child dictionaries, branches and node metadata remain. Continue using the recipient; the donor is neither an independently usable Tree nor a synchronized subtree view. Pass `donor.copy_treeobject()` to retain an independent original. An explicit connecting Branchstruct transfers directly too; copy it first when retaining an independent original. No automatic copy/reset, ownership guard, donor-cache synchronization or new rollback guarantee is added. Leaf strings create no separately maintained donor Tree.

A supplied `connect_branchstruct` defines the connection and takes precedence over `connect_length`; otherwise a fresh branch uses `connect_length`, defaulting to 0.0. A donor label, stem length or stem annotations never transfers implicitly to the connection or recipient root. To reuse a donor suffix deliberately, supply it as the connecting struct after setting an actual numeric structural length. Supplying `connect_length` too does not fill a missing length in that struct. The recipient retains its root/suffix, and donor root-node annotations stay with the surviving donor-root node. Detached basal branches remain separate and can also be supplied explicitly.

Successful grafts preserve the recipient's True/False/None interpretation at actual graft entry, for both Tree and leaf-name donors. Donor state, root degree, stem data and connecting length do not select a new interpretation. In a compound operation, graft preserves the state left by earlier steps; it does not restore an older historical state.

#### Complete SPR

Successful `spr()` preserves rootedness through its pruning/graft components. It retains the suffix if the old root survives and clears it if pruning removes that root. Moved internal edges and unrelated edges retain annotations; merged/split/connecting edges start with fresh non-length fields. Surviving node data follows final IDs, deleted data is removed, and the graftpoint starts empty. A moved leaf's node data is retained locally without changing extraction's leaf-string API. The connection uses detached basal length plus the other root-child length when pruning a direct child of a bifurcating root; its existing selection, placement and length calculations remain. The temporary donor is consumed. Preserving data does not recalculate its scientific meaning.

#### Direct numerical length edits

`setlength()` and native structural `set_branch_attribute(..., "length", value)` retain node/branch annotations and clear length-derived caches. Caller scaling loops have the same contract; statistical annotations are neither scaled nor recomputed automatically. `set_blens_from_heights()` reads explicitly stored node heights and writes represented lengths only. Rootedness, node IDs and the external root suffix remain unchanged. An extra annotation named `length` is separate data, not a structural write. Raw `.length` assignments require `clear_length_caches()` before derived queries; no watched fields or automatic annotation conversion is added.

`nodedist()` keeps computed distances on each Tree. Discarding a Tree releases its distance cache. Length or structural edits clear that Tree's distances, and `clear_length_caches()` does the same after raw length writes. Clearing one Tree does not clear another Tree's distances. Use `clear_length_caches()` in place of `nodedist.cache_clear()`. The method no longer provides the `functools` attributes `cache_info()` and `cache_parameters()`.

---

### 4. Building summary trees

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

# 3) Select output labels and metadata
print(sumtree.newick(
    labelfield="bipartition_cred",
    branch_attributes=("bipartition_cred", "length", "length_sd"),
    print_meta=True,
))
```

Parameters

- `treetype`: `"con"`, `"all"`, `"mcc"`, `"mbc"`, `"hip"`, `"mrhip"`
- `blen`: `"biplen"`, `"cladeheight"`, `"caheight"`, `"input"`, `"none"`
- `rooting`: `"mid"`, `"minvar"`, `"og"`, or `None`. None preserves the construction interpretation: consensus/all-compatible consensus/MBC are unrooted, while MCC/HIPSTR/mrHIPSTR are rooted. A successful explicit rooting operation establishes True.

Notes

- `"caheight"` needs access to the original tree sample again. Pass `count_burnin_filename_list`
  (the same structure used by `sumt`: tuples of `(count, burnin, filename)`), or do a second pass
  yourself with `CAHeightEstimator`.
- Select printing options when calling `Tree.newick()` or `Tree.nexus()` on the result.

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

### 5. Summary ownership

#### Independent input accumulation and constructed outputs

`TreeSummary.add_tree()` reads observed structural branch lengths into its own bipartition accumulators. It does not use input branches as accumulators or inherit sample labels, annotations or pre-existing statistics. Each newly observed bipartition gets one fresh record; later observations update that record using the existing numerical algorithms. Lazy summary finalization therefore leaves input branch data unchanged. Optional retained representative Trees have the separate storage contract below.

Consensus, compatible-bipartition insertion, MBC and branch annotation write-back materialize independently editable output branches and retained mutable annotations. Output branches contain the resulting lengths and intended annotations, with fresh accumulator/scratch fields and no quantile histograms. Repeated outputs are independent of one another and the summary. Annotation write-back can replace target Branchstruct objects; reacquire branch references from the returned Tree after a call. The low-level `from_biplist()` and `split_off_children()` direct-assignment contracts are unchanged. Exact-key interval write-back still preserves target entries that the calculation does not write.

#### Summary merging transfers ownership

Complete per-file statistics and diagnostics before `recipient.update(donor)`. Matching records combine into existing recipient records; unseen branch, clade, root and topology records may transfer directly, including histograms, interval dictionaries and child-clade-pair sets. Use the recipient afterward. Further donor accumulation, statistical reads or merging into another independent result are unsupported, because recipient updates/finalization can change transferred records. `len(donor)` remains available for immediate scalar progress reporting. No automatic donor copying, reset or reuse guard is added. This is separate from independent input accumulators and independently editable constructed outputs.

#### Retained representative Trees

With topology tracking and `store_trees=True`, the first Tree for each topology is retained by reference, potentially in both topology dictionaries. Leave its topology, names, lengths, rootedness and stored metadata unchanged while any summary retains it, including after a consuming merge. Read-only queries, printing and normal derived-cache work remain allowed. Edit an explicit `copy_treeobject()` instead. No automatic snapshot, mutation guard or ownership registry is added. The default `store_trees=False` imposes no representative-storage lifetime restriction. Representatives are retained input samples; their arbitrary annotations do not become aggregate statistics. Constructed summary outputs have a separate independent ownership contract.

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

MCC selects the observed rooted topology whose product of clade posterior frequencies is maximal, or equivalently whose sum of log clade frequencies is maximal. It reconstructs an independently editable tree from that clade topology and summary records; it does not return or edit a retained input representative.

This approach prioritizes topologies built from highly supported clades.

##### Maximum bipartition credibility (MBC)

The MBC tree is analogous but uses branch bipartitions instead of rooted clades. It therefore evaluates trees in terms of branch support rather than rooted clade support.

##### Practical difference

- MCC is usually the natural choice for rooted posterior trees, especially clock trees.
- MBC can be useful when input trees are unrooted or bipartition support is the quantity of interest.

Both methods choose an observed topology and construct an independent output. MBC uses an unrooted bipartition construction; MCC uses an oriented clade construction.

#### 9.3 HIPSTR and mrHIPSTR

`phylotreelib` also supports HIPSTR-style summary trees, following the idea of highest independent posterior subtree reconstruction.

The implementation recursively scores observed child-clade pairs and reconstructs an oriented clade hierarchy. HIPSTR can include clades below 50% support. mrHIPSTR favors clades with support above 50% and resolves the remaining hierarchy; it can also include lower-supported clades. Their rooted interpretation comes from clade construction, independently of root-frequency placement.

#### 9.4 Root summaries

When `trackroot=True`, `TreeSummary` can record where the root occurred across sampled trees.

`TreePostProcessor.root_maxfreq()` selects the most frequent recorded root bipartition compatible with the target tree, and a successful call establishes True. Its existing position calculation remains; unfinished length-fraction/rootcred work is separate. `set_rootcredibility()` annotates branches with `rootcred`, the posterior frequency with which each branch carried the root. Root summaries are useful when the posterior distribution contains uncertainty about root placement.

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

Ordinary Newick/Nexus labels belong to incoming branches. Reversing an existing edge retains its label and annotation values; separate merges and subdivisions follow the editing rules above. Directional labels and rooted-clade support can need scientific reassessment after rerooting.


---

## 10. Advanced notes: quantile tracking and the QuantileAccumulator

When `trackci=True`, `TreeSummary` tracks approximate quantiles for branch lengths and/or node heights. Also supply `ci_probs`, for example `(0.8, 0.95)`, to request central, equal-tailed credible intervals. With no probabilities specified, no intervals are requested. Enable `trackblen=True` for branch-length statistics or `trackheight=True` for clade-height statistics, alongside the corresponding `trackbips` or `trackclades` flag.

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

For distributions with meaningful variation, the method is approximate rather than exact.

The benefit is that it avoids storing all raw values while still providing useful credible intervals and quantiles for summary-tree annotation.

There is one exact special case for node-height summaries. `QuantileAccumulator` tracks the minimum and maximum observed values, including across merged accumulators. When `TreeSummary` or `CAHeightEstimator` determines that all observed heights are effectively identical, allowing for small numerical noise, it uses the observed height directly for the median and credible-interval endpoints instead of using a bucket midpoint. This keeps fixed node heights, such as fixed sampling dates at leaves, exact. Node-height distributions with real variation remain approximate, as do branch-length quantiles.

If exact quantiles are needed for a distribution with meaningful variation, they should be computed from the raw values outside the accumulator framework.

### 10.6 What the tracked quantiles are used for

Quantile summaries are typically used to populate metadata such as:

- branch-length credible intervals
- node-height credible intervals
- medians and related robust summaries

Print these with `Tree.newick()` or `Tree.nexus()`, passing `print_meta=True`. Exact annotation keys are needed only to restrict output. This complete example computes branch-length intervals for a consensus tree:

```python
import phylotreelib as pt

ts = pt.TreeSummary(
    trackbips=True, trackblen=True, trackci=True, ci_probs=(0.8, 0.95),
)
for sample in (
    "((A:1,B:1):2,C:3,D:3);",
    "((A:1.5,B:1.5):2.5,C:4,D:4);",
    "((A:1.2,B:1.2):2.8,C:4,D:4);",
):
    ts.add_tree(pt.Tree.from_string(sample))

sumtree = pt.build_sumtree(ts, treetype="con", blen="biplen")
sumtree.rooted = False  # No biological rooting was selected.
with open("consensus_intervals.nexus", "w") as out:
    out.write(sumtree.nexus(print_meta=True))

# Export only branch medians and 95% intervals, alongside ordinary branch lengths.
print(sumtree.nexus(
    print_meta=True,
    node_attributes=(),
    branch_attributes=("length_median", "length_95%_CI"),
))
```

The metadata names include `length_median`, `length_80%_CI`, and `length_95%_CI`. These branch-length quantiles remain approximate, including for constant lengths. For node-height intervals, see the [README example](../README.md#compute-credible-intervals-and-export-a-summary-tree) and the [common-ancestor-height recipe](recipes.md#build-an-mcc-tree-with-common-ancestor-heights-and-credible-intervals).

---

## 11. Public vs internal API

For everyday use, it is reasonable to treat the following as the stable public interface:

- `Tree`
- `Treefile`, `Newicktreefile`, `Nexustreefile`
- `TreeSet`
- `TreeSummary`
- `Distmatrix`
- `TreeError`
- `MISSING`
- `build_sumtree`

The remaining classes are useful, but many of them are better regarded as advanced helpers or implementation-level infrastructure unless they are explicitly needed by your workflow.

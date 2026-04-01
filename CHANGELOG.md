# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.2.1] - 2026 April 1

#### Fixed
- `Nodestruct` now includes the slot fields required by parsimony methods (`state`, `primary_set`, `secondary_set`, `optimal_set`, `fit`, `wasambig`), fixing failures when using `parsimony_possible_states`, `parsimony_assign_fits`, and `parsimony_count_changes` with slot-based node objects.

#### Changed
- Test suite reorganized and modernized: old tests migrated to pytest style, test files consolidated into a single combined test module, and test data moved into dedicated subdirectories under `tests/data/`.

---

## [2.2.0] - 2026 April 1

#### Added
- `Tree.randtree` now simulates proper Yule or coalescent tree topologies, producing more realistic tree structure. Default is to create an ultrametric tree with exactly the user-specified tree height. Setting `rate_sd > 0` adds lognormal noise to branch lengths, breaking ultrametricity.

#### Fixed
- NEXUS format output now consistently uses spaces for indentation (previously mixed spaces and tabs).

---

## [2.0.0] - 2026 March 17

This is a major release with breaking API changes. See the migration guide in `README.md` under "Upgrading from 1.x to 2.x" for full details.

#### Added
- Explicit summary-tree pipeline via `SummaryTreeBuilder` and `TreePostProcessor`, with optional `CAHeightEstimator` for common-ancestor heights and credible intervals.
- Credible interval support via `QuantileAccumulator` and the `trackci` / `ci_probs` parameters in `TreeSummary`.
- `PrintSpec` class plus helpers `configure_basic_printing(...)` and `configure_sumtree_printing(...)` for consistent, configurable Newick/NEXUS output.

#### Changed
- `Tree.graft(...)` now grafts onto a specified edge (parent → child) and uses `connect_length` / `connect_branchstruct` for the connecting branch.
- `Tree.spr(...)` updated to match the new grafting model and by default preserves total tree length.
- `Tree.subtree(...)` / `Tree.prune_subtree(...)` now return `(subtree, basal_branch)`; leaf subtrees return the leaf name as `str`.
- `Tree.reroot(...)` signature updated to use the same distance/fraction placement convention as other editing methods.
- `Tree.newick(...)` / `Tree.nexus(...)` now treat `None` as "use PrintSpec defaults".

#### Renamed (terminology: `depth` --> `height`)
All names relating to node heights — the distance from the tips to an internal node — have been renamed from `depth`-based to `height`-based, aligning with the convention in most phylogenetic software. **This is a breaking change; the old names are removed.**

| v1.x name | v2.x name |
|---|---|
| `TreeSummary(..., trackdepth=True, ...)` | `TreeSummary(..., trackheight=True, ...)` |
| `Tree.nodedepth(node)` | `Tree.nodeheight(node)` |
| `Tree.set_blens_from_depths()` | `Tree.set_blens_from_heights()` |
| `TreePostProcessor.set_mean_node_depths(tree)` | `TreePostProcessor.set_mean_node_heights(tree)` |
| `CADepthEstimator` | `CAHeightEstimator` |
| `build_sumtree(..., blen="meandepth", ...)` | `build_sumtree(..., blen="cladeheight", ...)` |
| `build_sumtree(..., blen="cadepth", ...)` | `build_sumtree(..., blen="caheight", ...)` |
| Node attributes `depth`, `depth_sd`, `depth_median` | `height`, `height_sd`, `height_median` |

#### Deprecated
- `Tree.deroot()` → use `Tree.collapse_bifurcating_root()`.

#### Fixed
- Multiple refactorings for improved performance (speed and memory) and bug fixes across tree editing, summary-tree construction, and output formatting.

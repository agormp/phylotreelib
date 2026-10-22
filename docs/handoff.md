# Development handoff

Updated: 2026-10-22. This document records project context, not permission to start deferred work.

## Read first

- The Nexus metacomment parser is planned but explicitly deferred by Anders. Do not implement it until asked. See [the parser plan](plans/nexus-metacomment-parser.md).
- Do not push, publish to PyPI, or change version numbers without a new request. `pyproject.toml` remains at 2.3.1 despite the unreleased changes in the checkout.
- The printing changes and initial handoff were fast-forwarded into local `main` at `15ed417`, then rootedness changes through `ef90599` were also merged. Both completed branches, `codex/explicit-tree-printing` and `codex/nexus-rootedness`, were deleted. The checkout is on local `main`; none of these changes has been pushed or released.
- Inspect `git status`, branches, and recent history at the start of a new session. Do not reset unrelated user changes or delete older experimental branches.

## Repositories

- phylotreelib: `/Users/agpe/Documents/code/python/phylotreelib`
- sumt: `/Users/agpe/Documents/code/python/sumt`
- Original investigation workspace: `/Users/agpe/Documents/4_archive/2026/2026_phylotreelib_debug`
- [sumt handoff](../../sumt/docs/handoff.md)

## Completed printing work

Relevant commits, in order:

- `24628c9`: removed stored printing configuration and made output arguments explicit.
- `c6a5e62`: documented the API break and migration.
- `69a6ab8`: automatic populated metadata export and seven-digit TreeSet precision.
- `86e50cc`: complete CI examples and documentation of the revised defaults.

Current contract:

- `Tree.newick()` and `Tree.nexus()` control output at the call site. There is no stored print specification.
- Removed APIs: `PrintSpec`, `configure_basic_printing`, `configure_sumtree_printing`, `Tree.set_print_spec`, and `Tree.get_print_spec`.
- `Tree` and `TreeSet` default to seven significant digits. Scalar `None` uses the fixed default.
- `print_meta=False` suppresses node/branch annotations. With `True`, an omitted or `None` selection includes all populated annotation fields; an empty sequence includes none; an explicit sequence includes only those fields. Node and branch selections are independent.
- Explicit attribute selections must include `"ci"` to emit intervals. `ci_labels=None` includes all stored intervals, `()` suppresses intervals, and explicit labels select those present.
- Automatic branch metadata excludes `mean`, `M2`, `n`, `quantiles`, `parent_height`, and `kid_height`. Zero and false values are retained. Empty values and missing fields are skipped.
- Strings are quoted, booleans use true/false, and collections use braces. Sets are sorted for repeatable output.
- Existing non-None struct defaults, such as `Branchstruct.length_sd = 0.0`, also count as populated; automatic output does not establish whether a statistic was estimated. Changing the data model for unset statistics is a separate decision.

README, CHANGELOG, docs/api.md, and docs/recipes.md cover the new behavior. Historical release notes intentionally retain descriptions of removed APIs as historical information.

## Credible intervals and previous investigation

- Both `trackci=True` and `ci_probs=(...)` are needed to request intervals. The old README omitted the probabilities; this is now fixed.
- Intervals are central/equal-tailed, not HPD intervals. Variable distributions use approximate log-bucket quantiles.
- The already-released constant-node-height fix preserves the observed height for median and CI endpoints and uses zero SD for effectively constant node-height distributions. Keep coverage of merged accumulators and genuinely variable leaf heights.
- Do not extend that special handling to `bipartsummary` without a separate design decision: non-positive branch lengths currently share a zero sentinel in quantile tracking.
- Rooting declarations, rerooting, and constant-height quantiles are distinct issues. `[&R]` communicates rootedness; it does not prevent a user or viewer from rerooting, nor does it recalibrate height annotations after rerooting. Do not infer that a particular user manually rerooted a tree merely from a mismatch.

## Verification baseline

The merged rootedness feature adds explicit `rooted=True/False/None` declarations to Tree/TreeSet Nexus output, with no inference from their convenience/traversal roots. See the [design decision](decisions/nexus-rootedness.md). After merging, the local main checkouts passed 270 phylotreelib tests and 123 sumt tests. Biopython 1.88 previously verified both emitted rootedness declarations. The parser is still deferred and versions are unchanged.

Eleven sumt CLI workflows were checked against the pre-rootedness main baseline; only the expected Nexus declaration changed. Serial and parallel output still matched, and Biopython recognized each Nexus output. The four CI documentation snippets were rerun successfully. See also the [dependency compatibility recommendation](../../sumt/docs/decisions/dependency-compatibility.md); no requirement bounds were changed yet.

Before the new rootedness work: 234 phylotreelib tests and 99 sumt tests pass. Four documentation snippets were executed, including clade-height, CA-height, and branch-length CI exports. The recipe's example download was replaced with the matching local test fixture during verification.

Four sumt CLI workflows (serial/parallel HIPSTR CA heights, MCC clade heights, consensus branch lengths) matched the pre-refactor library output exactly. Serial and parallel CA-height output also matched. sumt's printing tests passed against the pre-refactor library as well.

Reproduce the library tests from this repository:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD" python -m pytest -q -p no:cacheprovider
```

The interpreter used was Python 3.12.4 with pytest 8.3.2 and numpy 1.26.4. No installation is needed to test the sibling repositories; select them with PYTHONPATH. Do not rely on temporary verification scripts in `/private/tmp` surviving a restart.

## Working preferences and workflow

- Anders wants to understand consequential changes. Explain the approach before editing, keep changes scoped, and use meaningful conventional-style commit subjects.
- Use pytest, with one test class per method, in `tests/test_phylotreelib_combined.py` where appropriate.
- Prefer short-lived, descriptive branches from current `main`, e.g. `codex/nexus-rootedness` and later `codex/nexus-metacomment-parser`. Test, review, merge locally, then delete the completed branch. A permanent `dev` integration branch is unnecessary for one feature at a time.
- The same feature branch name can be used in both repositories for coordinated changes, but commits and merges remain independent.
- Branches such as `dateopt`, `constantheight`, and `parse_regex` predate this work. Do not delete or merge them without checking their contents and obtaining authorization. `dateopt` contains earlier date-option sketches, not a completed feature.
- A fresh chat is appropriate for the parser feature once authorized. Start with this document and the parser plan; the previous chat transcript is not required.

## Release checklist, not yet authorized

1. Complete and review the parser feature before the proposed phylotreelib 3.0.0 release. Anders wants a substantive new feature with the breaking print refactor.
2. Recheck migration examples, CI behavior, Nexus interoperability, and performance.
3. sumt now uses the new `rooted` export keyword. Update its minimum phylotreelib dependency when the release version is decided. Until then test against the matching local main checkout, not only the installed 2.3.1 wheel.
4. Choose versions and publish only after explicit approval. Current version numbers are not a declaration that unreleased code is compatible with the published packages.

# Maintenance follow-ups

Updated: 2026-10-28.

These items are deferred from the parser work. The original [parser plan](plans/nexus-metacomment-parser.md) and handoffs record the earlier decisions. The [API reference](api.md) describes current supported behavior.

## Caches and retained memory

- **Clustering scratch (C7):** `cluster_n()` stores a list that its next call recomputes. Review whether to keep that public scratch field or make it local.
- **Edge-count cache invalidation (C8):** length changes clear edge-count distances unnecessarily. This is excess invalidation, not a demonstrated stale-result bug.
- **Array caches (C9):** `pathdist_as_ndarray()` checks cache slots without populating them. Decide how cached mutable arrays should be owned before adding reuse or removing the slots.
- **Distance memoization (C11):** `nodedist()` uses a shared, unbounded function cache that can retain Trees and couples invalidation across instances. Measure representative workloads and lifetimes before choosing a bounded or per-instance replacement.
- Review the vestigial `_frozenset_leaves`/`_leaf2index` fields and unused `clear_caches(preserve=...)` interface separately from correctness fixes.

## Tree operations

- **Subtree iteration:** `Tree.__iter__()` visits leaves and tries to attach `.basalbranch` to the leaf-name string returned by extraction. Repair needs an explicit decision about yielded values. Keep the existing leaf-extraction interface until that separate contract is agreed. This concerns iteration over subtrees of one Tree, not iteration over a treefile.
- **Root-frequency selection:** the current API/state repairs do not complete `root_maxfreq()`'s documented length-fraction adjustment or root-credibility assignment. Review intended numerical behavior and verification separately.

## Input scope and future workflows

- Broader quoted taxon/internal/Nexus names, underscore/case changes, additional dialects and tree weights remain outside the supported parser scope.
- Summarizing caller-chosen imported annotations needs its own design. Current summary observations come from topology and structural lengths; imported sample annotations are not automatically aggregate statistics.
- TreeAnnotator fixture provenance, reruns and possible upstream reports remain separate investigations. Existing parser fixtures do not establish the scientific correctness of producer statistics.
- The separate `--date` feature remains future work.

## Numerical reproducibility

- Serial and parallel summaries have differed by about `4.1e-16` in node-height SDs around `2.47e-13`, including in the pre-parser code. Some repeated parallel runs also differ in their final printed digits. Investigate stable accumulation and merge order with the existing comparison rule unchanged.

## Before release

Release preparation targets phylotreelib 3.0.0 and sumt 4.2.0, with Python 3.10 or newer and `phylotreelib>=3.0.0,<4`. Verify the built distributions together in fresh environments before publication. Source-checkout tests do not establish compatibility of published packages. See sumt's dependency-compatibility decision for the release sequence.

# Maintenance follow-ups

Updated: 2026-10-28.

These items are deferred from the parser work. The [historical parser plan](https://github.com/agormp/phylotreelib/blob/4ffca97637156a992b8c80c191ce0fd763a21e3e/docs/plans/nexus-metacomment-parser.md) records the initial design. The [API reference](api.md) describes current supported behavior.

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

## Numerical precision

Version 3.0.0 retains the rounding remainder when updating online means. This corrects the near-zero node-height SD discrepancy found between serial and parallel summaries.

## Before release

Release preparation targets phylotreelib 3.0.0 and sumt 4.2.0, with Python 3.10 or newer and `phylotreelib>=3.0.0,<4`. Verify the built distributions together in fresh environments before publication. Source-checkout tests do not establish compatibility of published packages. See sumt's dependency-compatibility decision for the release sequence.

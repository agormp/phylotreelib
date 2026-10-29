# Maintenance follow-ups

Updated: 2026-10-29.

This lists remaining follow-ups and work completed since the parser release. The [API reference](api.md) describes current supported behavior.

## Caches and retained memory

- **Clustering scratch (C7):** `cluster_n()` stores a list that its next call recomputes. Review whether to keep that public scratch field or make it local.

## Tree operations

- **Root-frequency selection:** the directly called `root_maxfreq()` helper still does not complete its documented length-fraction adjustment or root-credibility assignment. Review its intended behavior separately. This is not a demonstrated defect in the published HIP/CA workflow.

## Input scope and future workflows

- Broader quoted taxon/internal/Nexus names, underscore/case changes, additional dialects and tree weights remain outside the supported parser scope.
- Summarizing caller-chosen imported annotations needs its own design. Current summary observations come from topology and structural lengths; imported sample annotations are not automatically aggregate statistics.
- TreeAnnotator fixture provenance, reruns and possible upstream reports remain separate investigations. Existing parser fixtures do not establish the scientific correctness of producer statistics.
- The separate `--date` feature remains future work.
- Using medians as sumt point estimates remains future work. Median annotations are already available with `--ci`.

## Completed since the 3.0.0 release (unreleased)

- Tree iteration yields proper internal subtrees with copied basal branches. The source root and individual leaves are excluded.
- Node-distance caches belong to each Tree. Discarded Trees can be collected, and clearing one Tree leaves other Trees' distances intact.
- Internal-node views are no longer stored in Trees, so queries do not prevent pickling or deep copying.
- Length-only edits retain cached edge counts. Structural changes and renaming clear them. Arrays remain fresh and writable, and unused private cache fields have been removed.
- `clear_caches(preserve=...)` retains its public signature. Listed structural caches are preserved, and length caches are always cleared.

## Numerical precision

Online mean and SD calculations use floating-point approximations. Serial and parallel SDs can differ near zero. Comparisons allow for rounding at the scale of the underlying heights or branch lengths.

## Release status

phylotreelib 3.0.0 and sumt 4.2.0 were published and verified together on 2026-10-28. They require Python 3.10 or newer, and sumt requires `phylotreelib>=3.0.0,<4`. The completed work listed above is not yet released. Verify built distributions together in fresh environments before a future publication. See sumt's dependency-compatibility decision for the release sequence.

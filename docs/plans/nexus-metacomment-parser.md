# Nexus metacomment parser plan

Status: agreed direction, implementation explicitly deferred on 2026-10-22.

## Goal and scope

Read typed Nexus-style tree annotations, including annotations on the root, and export them again without losing unknown keys. Reuse the topology state machine in `NewickStringParser` through a derived annotated-Newick parser. Preserve a fast path for ordinary Newick.

This is the proposed substantive feature for phylotreelib 3.0.0, accompanying the breaking printing simplification. It is not permission to implement the parser or bump versions. The separately requested rootedness-output feature does not constitute implementation of this parser.

## Findings in the current code

- `NewickStringParser.__init__` already calls overridable `set_token_delimiters()` and `update_dispatch()` hooks. Dispatch entries contain a handler and the next state.
- `parse()` removes all whitespace and then splits on delimiters with a regular expression. This destroys quoted whitespace and splits punctuation within annotation values. Merely adding `[` and `]` to delimiters is insufficient.
- `Tree.from_string()` selects the Newick parser, while `Tree._from_string_private()` already supports reuse of a parser instance by file readers.
- `TreefileBase.get_treestring()` splits on raw semicolons and reconstructs the remaining buffer without them. It also moves everything after the last `)` to `below_root`, removing root labels/metadata before parsing.
- `Nexustreefile.__init__()` strips comments while reading the first tree's header, then discards the prefix before `(`, including rootedness flags. Subsequent trees follow a different path. Both paths must preserve annotations consistently.
- Translate-block parsing splits on commas and whitespace, so quoted names require attention too.
- `Nodestruct` and `Branchstruct` use `__slots__`; arbitrary imported keys cannot simply be assigned with `setattr`.
- All `Tree` objects have an implementation root, including representations of unrooted trees. Root degree is not a reliable declaration of biological rootedness.

## Proposed implementation batches

### 1. Tokenization extension

Add a tokenizer hook to the base parser while retaining its ordinary-Newick behavior. The derived parser should scan quoted strings, balanced comments, and brace-delimited collections, yielding complete comment tokens alongside topology tokens. Preserve whitespace inside values and report malformed input with useful locations. Do not use `eval`.

Keep tokenization separate from topology construction. Comments should generally be self-loop transitions in the appropriate state, not a separate replacement topology parser.

### 2. Typed metadata grammar and storage

Support integers, signed/scientific floating-point values, booleans, quoted and simple unquoted strings, and brace-delimited collections. Distinguish ordinary comments from `[&key=value,...]` metacomments. Specify quote escaping and missing-value behavior against real producer fixtures before coding.

Add per-node/per-branch metadata mappings for unknown keys and a tree-level location for tree annotations. Decide the authoritative storage and lookup rules for recognized fields so that two copies cannot diverge. Existing `height`, median, SD, and `.ci` fields should remain usable through the current API.

Preserve unknown annotation names and types. Do not rename HPD intervals to CI intervals, or silently overwrite the actual branch length with a comment's `length` value. Imported height annotations must remain distinct from the branch-derived `nodeheightdict` cache.

### 3. State-machine extension

Add comment transitions after leaf names, internal-node closures, labels, and around branch lengths. Handle trailing root annotations without assuming an incoming root branch exists. Define root-label and stem-length storage explicitly, preserving any required `below_root` compatibility.

Define node-versus-branch ownership as a documented import/export convention. Our writer places node annotations before the colon and branch annotations after the length; other producers can group differently named statistics in one comment. Do not infer all ownership from key names. Also consider ambiguous placement when branch lengths are suppressed on export.

Handle `[&R]` and `[&U]` as tree-level rootedness declarations, not root-node attributes. Later importing rootedness must not reintroduce hidden formatting defaults: export arguments remain authoritative.

### 4. File-reader integration

Use quote/comment-aware statement termination rather than raw semicolon splitting. Preserve annotation text through header processing and root handling. Treat the first and subsequent trees identically, including multiple tree statements on one line. Make translation tables quote-aware.

Proposed entry points: `Tree.from_string(..., parse_meta=True)`, `Treefile(..., parse_meta=True)`, and corresponding concrete readers. Keep metadata parsing opt-in initially. Define interaction with `strip_comments` so requesting parsed metadata does not simultaneously discard it; ordinary comments can still be ignored.

### 5. Printing and interoperability

Extend automatic and selective metadata output to imported mappings using the same `node_attributes`, `branch_attributes`, `ci_labels`, and `print_meta` contract. Ensure escaping is compatible with the new parser and verify representative external consumers. Round-trip semantic values; exact whitespace reproduction need not be promised.

Initial scope should be Nexus/BEAST-style metacomments, not every possible extended-Newick dialect. Explicitly document whether NHX and other producer-specific conventions are unsupported.

### 6. Tests, benchmarks, and sumt integration

Use one pytest class per method, extending the existing combined module. Cover:

- Scalar types, collections, nested values, escaped quotes, whitespace, punctuation, and missing values.
- Leaf/internal/root annotations, annotations with/without lengths and labels, repeated comments, unknown keys, and CI versus HPD names.
- Rooted/unrooted declarations, translated/quoted taxa, first/subsequent trees, multiline comments, multiple trees per line, and semicolons in quoted values.
- Malformed annotations, duplicate-key policy, accurate error context, and safe handling of arbitrary text.
- Reusing parser objects without leaking annotations between trees.
- Copying and tree-editing operations, so metadata mappings are not accidentally shared or lost.
- Semantic parse/export/parse round trips of phylotreelib output and real BEAST/sumt fixtures.
- Existing plain-Newick regression tests and a benchmark of large posterior samples with metadata parsing disabled and enabled.

sumt should continue deriving summary heights from branch lengths. It should not copy arbitrary annotations from one posterior sample onto the summary tree. Enable metadata input only for a workflow that actually needs it. A future `--date` feature remains separate and was not implemented as part of printing work.

## Decisions to resolve before coding

1. Authoritative storage and collision policy for recognized and unknown keys.
2. Node/branch annotation placement conventions, including exports without branch lengths.
3. Duplicate-key handling, missing-value representation, and supported quote escaping.
4. Root label/stem length/`below_root` compatibility and explicit imported rootedness storage.
5. Precise interaction between `parse_meta` and `strip_comments`.

## Suggested workflow

Start from tested local `main` only after authorization, on `codex/nexus-metacomment-parser`. Split commits into tokenizer hooks, typed parsing/storage, reader integration, round-trip printing, and documentation/tests as appropriate. Preserve the printing simplification. Do not publish or change versions until the release has been approved.

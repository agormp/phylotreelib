# Explicit Nexus rootedness declarations

Date: 2026-10-22. Status: implemented, tested, and merged into local `main` in both repositories. The completed `codex/nexus-rootedness` branches were deleted. Not pushed or released.

## Context

Every phylotreelib Tree is internally rooted for traversal. When parsing Newick, the outermost parentheses identify that implementation root. As Anders clarified, this is a convenience rather than a statement of biological rootedness. An unrooted topology still needs a starting node for this representation. A rooted tree can also have more than two root children, so root degree is not a safe test.

The existing printing refactor deliberately removed stored formatting configuration. Rootedness output should remain explicit at the export call, not become another hidden setting or an inference from traversal state.

## Decision

Append `rooted=None` to `Tree.nexus()` and `TreeSet.nexus()` without moving existing positional arguments:

- `True`: emit `[&R]` immediately before the Newick tree description.
- `False`: emit `[&U]` at the same position.
- `None`: omit the declaration, preserving the previous default output.
- Reject other values, including integer 0/1 and strings, with TreeError.

Example: `tree nexus_tree = [&R] (A:1,B:1,C:1);`

The declaration is independent of `print_meta`. Suppressing node/branch statistics must not silently change the declared interpretation of a tree. It never moves the root, changes branch lengths, alters node annotations, or persists a printing setting. TreeSet applies the one supplied declaration to all its trees. Newick methods are unchanged.

Omitting the declaration is not an explicit statement that a tree is unrooted. The Nexus specification describes rooted as the default; consumers may have their own heuristics, making an explicit flag preferable when the caller knows the intended interpretation.

The current input readers accept these comments by discarding them, but do not retain rootedness information. Importing and round-tripping that information belongs to the deferred metacomment-parser project.

## sumt integration

`print_sumtree()` has workflow information that a generic Tree lacks:

- MCC, HIPSTR, and mrHIPSTR output is declared rooted.
- Consensus, all-compatible consensus, and MBC output is declared rooted only when an explicit rooting workflow was chosen; otherwise it is declared unrooted.
- `args.actively_rooted` can be an outgroup list rather than a bool; normalize it before passing the export argument.
- Do not infer rootedness merely from a height-estimation option.
- Keep the declaration with `--nometa`; do not add it to plain Newick output.

This does not stop a downstream user from rerooting or changing a viewer's settings. It does not itself resolve stale annotations after rerooting or quantile approximation issues.

## Compatibility evidence

Primary sources accessed 2026-10-22:

- [Maddison et al., Nexus format specification, TREES block](https://www.phylo.org/sub_sections/nexus_document.pdf): specifies `[&R]` and `[&U]` at the beginning of the tree description. This is a standard declaration, not a sumt-specific annotation.
- [MrBayes 3.2 manual, section 5.2](https://raw.githubusercontent.com/NBISweden/MrBayes/develop/doc/manual/Manual_MrBayes_v3.2.pdf): documents these rooted/unrooted declarations. This does not imply support for arbitrary node metadata; the manual separately limits such input.
- [treeio's BEAST example](https://github.com/YuLab-SMU/treeio/blob/devel/inst/extdata/BEAST/beast_mcc.tree): contains a leading `[&R]` alongside BEAST node metadata, showing the convention in an established interoperability fixture.

Local check: installed Biopython 1.88 reads emitted `[&R]`/`[&U]` and reports the corresponding rooted boolean. The phylotreelib tests check parseability, translated names, multiple trees, metadata suppression, no inference from root degree, and call independence.

Inference: the declaration is appropriate for standard Nexus consumers. No claim is made that every viewer/version was tested. FigTree and MrBayes graphical/runtime behavior have not been exercised here; a manual FigTree smoke test before release remains useful.

## Alternatives considered

- Always emit `[&R]`: rejected because an implementation root does not establish biological rootedness.
- Infer rootedness from the number of root children: rejected because rooted polytomies and arbitrary traversal roots make this ambiguous.
- Store a new inferred boolean on every Tree now: deferred until input semantics and metadata ownership are designed with the parser feature.
- Gate the declaration on `print_meta`: rejected because rootedness is interpretive information, not optional statistical annotation.

## Release gate

Version numbers remain unchanged. This sumt feature requires a phylotreelib version that accepts the new keyword. Test the sibling source checkouts together, and raise sumt's dependency only when the compatible release number is chosen. Do not publish with the unchanged minimum dependency. The full metacomment parser remains deferred.

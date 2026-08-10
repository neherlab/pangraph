# Graph merging — design manifesto

Status: **`pangraph merge` implemented; verification and `build`-side name checks still pending (see §6)**
Integration branch: `feat/merge` (merged into `master` last, after all phase branches below)

This document describes the design for a new `pangraph merge` command, which combines two
pre-existing pangenome graphs into one. It covers the motivation, the identifier model, the
required changes to existing code, and the phased roadmap.

An accompanying `dev/design/merge-implementation-notes.md` will record as-built decisions and
gotchas as the phases land; this file describes the *intent*.

---

## 1. Motivation

Today `pangraph build` can only construct a graph from scratch, out of a set of FASTA sequences.
There is no way to add new genomes to an existing graph: the only option is to rebuild everything
from the full sequence set.

The proposed workflow for extending a graph is a two-step composition of existing pieces:

```bash
# 1. turn the additional genomes into their own (possibly singleton) graph
pangraph build new_genomes.fa -o new.json

# 2. merge it into the existing graph
pangraph merge base.json new.json -o extended.json
```

Step 2 is exactly the operation that already happens at every internal node of the guide tree
during `build`. The command is a thin CLI wrapper around it.

### What this buys us

- Adding genomes to a large graph without redoing the whole build.
- Merging graphs built from different sources or with different parameters.
- A directly testable, user-visible entry point into the merge machinery, which is currently only
  reachable through the full build pipeline (or the undocumented `bin/merge_two_graphs.rs`).

### What it does not buy us

`merge(build(A), build(B))` is **not** expected to be identical to `build(A ∪ B)`. The guide tree
and therefore the order in which homologous regions are discovered differ, so the resulting block
decompositions will differ in detail. Both are valid pangenome graphs over the same sequences.

---

## 2. Current state

### 2.1 What already exists

The core operation is complete and generic. `merge_graphs(left, right, args)`
(`packages/pangraph/src/pangraph/graph_merging.rs:26`) joins two graphs, iteratively self-merges
homologous blocks until convergence, and removes transitive edges. It makes **no assumption** that
its inputs are singletons or that they come from the same build. `bin/merge_two_graphs.rs` is a
development binary that already tries to call it directly.

Distance computation also already works on arbitrary graphs: `mash_distance` takes `&[Pangraph]`
and sketches all block consensuses (`distance/mash/mash_distance.rs:9`).

So this is primarily a plumbing and identifier problem, not an algorithmic one.

### 2.2 What is missing

1. **Identifier collisions** between independently built graphs (§3). This is the hard blocker.
2. `build` **panics on a single input sequence**, so the "wrap my new genome in a singleton graph"
   half of the workflow does not work (§4.1).
3. The merge parameters are entangled with `PangraphBuildArgs`, which carries build-only fields
   (§4.2).
4. **Path names are not guaranteed unique**, although the rest of the codebase treats them as the
   de-facto identity of a genome (§4.3).
5. Sequence verification is keyed on input order / FASTA index, which is meaningless for a merged
   graph (§4.4).
6. There is no `merge` subcommand (§4.5).

---

## 3. The identifier model

### 3.1 The problem

`Pangraph::singleton` (`pangraph/pangraph.rs:29`) assigns

```rust
let node_id  = NodeId(fasta.index);
let block_id = BlockId(fasta.index);
let path_id  = PathId(fasta.index);
```

i.e. small sequential integers taken from the position of the record in the input FASTA. Every
later identifier is content-derived, **but only for entities that are actually touched by a
merge**. Blocks that never find a homologous partner keep their original small integer id all the
way to the final graph.

Measured on two graphs built independently from disjoint subsets of `data/ges-1.fa` (6 + 6
sequences):

```
paths   A: 6    B: 6    overlap: 6   → {0,1,2,3,4,5}
blocks  A: 43   B: 33   overlap: 0
nodes   A: 94   B: 85   overlap: 0
```

and on a graph built from two unrelated sequences (influenza + mpox, no homology found at all):

```
blocks: {0, 1}    nodes: {0, 1}    paths: {0, 1}
```

So: **path ids always collide**, and **block and node ids collide whenever a block never merged** —
which covers every singleton graph, every accessory block, and every divergent genome.
`graph_join` (`graph_merging.rs:74`) panics on any duplicate key, so `merge` fails immediately
today.

### 3.2 Two properties of the current model

**Identifiers are labels, not content hashes.** `circularize::merge_blocks` deliberately *reuses*
an existing `BlockId` for a newly concatenated block (`circularize/merge_blocks.rs:110-142`).
Nothing in the codebase relies on a block id being a hash of its content, so we are free to
relabel entities as long as all references are rewritten consistently.

**Identifiers are re-hashed during the merge**, and those hashes consume the ids we relabel:

| Site | New id derived from |
|---|---|
| `pangraph_node.rs:45` | `(block_id, path_id, strand, position)` |
| `reweave.rs:132` | `(qry.name, qry.interval, reff.name, reff.interval)` |
| `pangraph_interval.rs:129` | `(block_id, interval)` |
| `detach_unaligned.rs:99` | `(node_id, seq)` |
| `pangraph_path.rs:44` | `(nodes, tot_len, circular, desc, name)` |

This is safe **provided the inputs are made disjoint first**. It also means a half-measure — for
instance deduplicating keys only inside `graph_join` — would not work: colliding inputs would
produce colliding *derived* ids later in the same merge.

### 3.3 `PathId` is semantically load-bearing

Unlike the other two, `PathId` currently carries meaning beyond identity: `reconstruct` sorts paths
by id and uses `path_id.0` as the output `FastaRecord.index`
(`commands/reconstruct/reconstruct_run.rs:56-76`), which `build --verify` then uses to index back
into the input FASTA records (`commands/build/build_run.rs:44`). For a graph produced by `build`,
`PathId` == input FASTA index, so reconstruction order == input order.

That correspondence cannot survive a merge (§4.4), but path ids should still be small, contiguous
and ordered rather than becoming arbitrary hashes.

### 3.4 Decision: relabel everything after the first graph

Before joining, the right-hand graph is relabeled so that its identifier space is disjoint from
the left-hand one. The **left graph is never modified**.

```rust
/// Re-derives every block, node and path id of the graph. Block and node ids become
/// `id((salt, old_id))`; path ids are renumbered contiguously from `path_id_offset`,
/// preserving their relative order. Consensuses, edits, names, descriptions, strands
/// and positions are unchanged.
fn relabel_in_place(&mut self, salt: usize, path_id_offset: usize)

/// Makes `right` id-disjoint from `left`. `left` is left untouched.
pub fn make_disjoint(left: &Pangraph, right: &mut Pangraph) -> Result<(), Report>
```

Rationale for relabeling *all* of the right graph's block and node ids, rather than only the ones
that actually clash: it is uniform, trivial to reason about, and trivial to test. The alternative
(remap-on-collision) would preserve more identifiers across an append, but block ids of an appended
graph are of limited interest to users compared with those of the base graph — and those are
preserved either way.

**Consequences, to be documented in the user docs:**

- Put the **base graph first**. Its block, node and path ids survive the merge unchanged (except
  where blocks are genuinely reweaved); the second graph's ids always change.
- `merge` is therefore *not* symmetric in identifiers, even though it is essentially symmetric in
  content.

### 3.5 Implementation notes

- **Relabel in place.** Block consensuses are the bulk of a graph's memory; building a relabeled
  copy would double peak usage on a large base graph. Use `std::mem::take` on the three maps and
  rebuild them in place.
- No new constructors are needed: `PangraphBlock::new`, `PangraphNode::new(Some(id), …)` and
  `PangraphPath::new(Some(id), …)` all accept an explicit id.
- References to rewrite: the three map keys, `PangraphBlock::alignments` keys,
  `PangraphNode::{block_id, path_id}`, and `PangraphPath::nodes`.
- **Assert injectivity.** After relabeling, the three maps must have the same lengths as before (a
  hash collision inside one graph would silently drop entities) and the resulting id sets must be
  disjoint from the left graph's. `make_disjoint` bumps the salt and retries if not; this is
  astronomically unlikely with `XxHash64`, but the check is cheap.
- Path ids are assigned as `path_id_offset + rank`, where `rank` is the position in the graph's
  existing (sorted) path key order, and `path_id_offset` is `max(left path ids) + 1`. This keeps
  ids contiguous and keeps the right graph's genomes in their original relative order, after all of
  the left graph's.
- Relabeling is deterministic: `utils::id::id` uses `XxHash64` with a fixed seed.

---

## 4. Required changes to existing code

### 4.1 `build` must accept a single sequence

```
$ pangraph build one.fa -o one.json
The application panicked: index out of bounds: the len is 1 but the index is 1
Location: packages/pangraph/src/tree/neighbor_joining.rs:28
```

`build_tree_using_neighbor_joining` (`tree/neighbor_joining.rs:16-33`) skips the
`while nodes.len() > 2` loop for a single graph and then indexes `nodes[1]`. An empty input panics
on `nodes[0]`.

Fix: return the lone leaf clade when `nodes.len() == 1`, and a clean error for an empty input. The
rest of the pipeline already tolerates a leaf-only tree — `postorder` no-ops on leaves, the root's
`data.take()` returns the graph, and `ProgressBar::new` already guards `n_total <= 1`. Worth a
companion test for `build_tree_from_newick` on a single-leaf Newick string, which shares the shape.

This is a prerequisite for the whole workflow: appending a single genome requires building a
one-sequence graph.

### 4.2 Extract `GraphMergeParams` from `PangraphBuildArgs`

`merge_graphs` takes `&PangraphBuildArgs`, and that type is threaded all the way down into
`MergePromise::solve_promise` (`reweave.rs:40`), `reconsensus_graph` (`reconsensus/reconsensus.rs:35`),
`edit_consensus_and_realign` (`pangraph_block.rs:295`) and `map_variations`
(`align/map_variations.rs:43`). Only five fields are ever used; the rest (`input_fastas`,
`output_json`, `circular`, `verify`, `guide_tree`, `no_progress_bar`) are build-only and meaningless
for `merge`.

```rust
pub struct GraphMergeParams {
  pub aln_args: AlignmentArgs,            // #[clap(flatten)]
  pub alignment_kernel: AlignmentBackend,
  pub max_self_map: usize,                // default 100
  pub extra_band_width: usize,            // default 5
  pub max_alignment_attempts: usize,      // default 4
}
```

`PangraphBuildArgs` and `PangraphMergeArgs` both `#[clap(flatten)]` it, so the two commands expose
an identical set of alignment options.

Place it next to `AlignmentArgs` in `align/`, and move `AlignmentBackend` there as well. Side
benefit: this removes the current `align → commands::build` dependency
(`map_variations.rs:3`, `reconsensus.rs:1`, `pangraph_block.rs:2`), which is backwards.

Mechanical follow-ups: the clap `default_value_t = PangraphBuildArgs::default().alignment_kernel`
becomes `GraphMergeParams::default()…`, and tests constructing `PangraphBuildArgs::default()` purely
to reach alignment parameters switch to `GraphMergeParams::default()`. Also generalize
`build_cmd_preliminary_checks` (the "is `mmseqs` on PATH" check) to take `&GraphMergeParams`, so
`merge` inherits it.

### 4.3 Path names must be unique

Path names are already treated as the identity of a genome throughout the codebase:
`Pangraph::path_id_by_name` (`pangraph.rs:258`), `simplify` (which unwraps `path.name()` outright,
`simplify_run.rs:26`), and the export commands. Nothing enforces uniqueness, however —
`build_run.rs:72` still carries a `// TODO: check for duplicate fasta names`. The only place that
errors today is `build_tree_from_newick` (`tree/newick.rs:88`), i.e. only when `--guide-tree` is
used.

Decision: **duplicate path names are a hard error, in both `build` and `merge`.**

- `build`: after reading the input FASTA records, error if two records share a name. This makes the
  existing `--guide-tree` behaviour unconditional and fills the TODO.
- `merge`: error if a name occurs more than once within either input graph or across the two.
  This is what catches the common mistakes of merging a graph with itself, or re-adding a genome
  that is already present in the base graph — which would otherwise silently produce two
  identically named paths and break `export` and `simplify` downstream.

Error messages should list the offending names, not just the first one.

Paths with no name at all are left alone by the duplicate check (multiple unnamed paths are not an
error), but `--verify` requires every path to be named, because names are the verification key
(§4.4).

### 4.4 Reconstruction and verification must be name-keyed

`reconstruct` emits one FASTA record per path, sorted by `PathId`, with `index = path_id.0`
(`reconstruct_run.rs:56-76`). For a graph produced by `build`, `PathId` is the input FASTA index, so
the output order reproduces the input order exactly and `index` is meaningful.

**Neither property can hold for a merged graph.** The right graph's paths are renumbered (§3.4), so
while the output stays deterministic (still sorted by path id), its order corresponds to no
original input FASTA ordering, and `index` is merely a position within the merged graph.
Reconstruction itself remains fully supported — every sequence comes back byte-identical — but
*order is not guaranteed*, and neither order nor `index` may be used to pair records up.

This has three consequences:

1. **`compare_sequences` must stop comparing whole records.** It currently tests
   `left != right` on `FastaRecord` (`reconstruct_run.rs:45-54`), and `FastaRecord` derives
   `PartialEq` over all fields *including `index`* (`io/fasta.rs:17-24`) — despite an error message
   that only mentions length. It should compare sequence contents and report the path name.

2. **Verification is keyed by path name, not by index or position.**

   ```rust
   fn verify_graph_sequences(graph: &Pangraph, expected: &BTreeMap<String, Seq>) -> Result<(), Report>
   ```

   `build` fills the map from the input FASTA records; `merge` fills it by reconstructing each
   input graph *before* merging. Both are now sound because names are unique (§4.3). The map form
   also keeps working for the intermediate graphs checked inside the build loop, which contain only
   a subset of the paths.

3. **`pangraph reconstruct` documentation must state** that record order matches the original input
   FASTA order only for graphs produced directly by `build`, and that consumers should match records
   by name.

### 4.5 The `merge` command

New module `commands/merge/{mod.rs, merge_args.rs, merge_run.rs}`, dispatched from
`commands/root_args.rs:61` and `commands/main.rs:19`.

```
pangraph merge <LEFT_GRAPH> <RIGHT_GRAPH> -o <OUTPUT>
```

`PangraphMergeArgs`: two positional graph paths, `-o/--output-json` (default `-`),
`-f/--verify`, and a flattened `GraphMergeParams`.

`merge_run` flow:

1. Load both graphs (`Pangraph::from_path`); `sanity_check` each under `debug_assertions`.
2. Preliminary checks: both graphs non-empty; the `mmseqs`-on-PATH check when that kernel is
   selected; duplicate path names within and across the inputs (§4.3); a **warning** if the two
   graphs disagree on circularity (`circular` is per-path, so mixing is structurally fine, but
   `build -c` is global and easy to get wrong for the new-genomes graph).
3. `make_disjoint(&left, &mut right)` (§3.4).
4. If `--verify`: reconstruct the expected sequences from both inputs into a
   `BTreeMap<String, Seq>` keyed by path name (§4.4).
5. `merge_graphs(&left, &right, &args.merge_params)`.
6. If `--verify`: `verify_graph_sequences(&merged, &expected)`.
7. Write the output JSON.

`bin/merge_two_graphs.rs` is superseded by this command — it cannot work today anyway, since it
calls `graph_join` on colliding ids — and should be removed in the same phase.

Exactly two input graphs are supported. See §7 for why N-way is deferred.

---

## 5. Cost model

`self_merge` runs an all-vs-all alignment over *every* block consensus in the joined graph on each
round. Appending 5 genomes to a 500-genome graph therefore costs roughly what the root merge of a
505-genome build costs: it is **not** proportional to the number of new genomes.

This is inherent to the existing algorithm and is no worse than what `build` already does at the
root of the guide tree, but users will reasonably assume that appending is cheap. It must be stated
explicitly in the user documentation.

---

## 6. Phased roadmap

Each phase is its own branch and PR, targeting the integration branch `feat/merge`. `feat/merge` is
merged into `master` last, once all phases have landed.

| # | Branch | Scope | Status |
|---|---|---|---|
| 1 | `feat/merge-single-seq-build` | §4.1 — NJ tree for 1 and 0 input graphs | landed |
| 2 | `feat/merge-cmd` | §4.2 — extract `GraphMergeParams`, no behaviour change | landed |
| 3 | `feat/merge-cmd` | §3.4 — `relabel_in_place` / `make_disjoint_from` | landed |
| 4 | `feat/merge-cmd` | §4.5 — the command itself, plus integration tests | landed |
| 5 | `feat/merge-unique-names` | §4.3 — duplicate path names are an error in `build` too | todo |
| 6 | `feat/merge-verify` | §4.4 — name-keyed verification shared with `build` | todo |
| 7 | `feat/merge-docs` | §9 — tutorial, `reconstruct` docs, CHANGELOG | todo |

Phases 2–4 were implemented together, since a `merge` command without §3.4 panics on the first
identifier collision and would not be testable.

Two pieces of §4.3 and §4.4 are therefore only half-done, and phases 5 and 6 finish them:

- Duplicate path names are rejected by `merge` (across and within its two inputs), but `build` still
  accepts duplicate FASTA names except when `--guide-tree` is used.
- Name-keyed verification exists as `verify_merged_sequences` inside `merge_run`, private to the
  merge command. `build` still verifies against input FASTA records by index, and
  `compare_sequences` still compares whole `FastaRecord`s (including `index`). Phase 6 unifies the
  two behind a single `BTreeMap<String, Seq>`-based helper.

---

## 7. Non-goals (for now)

- **N-way merge.** `merge` takes exactly two graphs. Merging N graphs along a guide tree is a
  natural generalisation — `mash_distance` already works on arbitrary graphs, so the neighbour-joining
  machinery would work unchanged — but it would require factoring the tree-traversal loop out of
  `build()`, and users can chain pairwise merges in the meantime.
- **FASTA input to `merge`.** Once `build` accepts a single sequence (§4.1), the two-step workflow
  composes cleanly. Accepting FASTA directly would duplicate `build`'s reading, `--circular`
  handling and name checks inside `merge`.
- **Guide-tree-driven merge order.** Only relevant for N-way.
- **Incremental performance.** See §5; making appends proportional to the number of new genomes
  would require a different alignment strategy and is out of scope.
- **Provenance metadata.** The graph JSON has no place to record which inputs a merged graph came
  from. Not added here.

---

## 8. Testing strategy

**Unit**

- `relabel_in_place` / `make_disjoint`: two hand-built graphs with identical ids `0..2` become
  joinable, and `sanity_check` passes on the join.
- Relabeling is an isomorphism: sequences reconstructed from the graph, keyed by path name, are
  unchanged by relabeling.
- Path ids after relabeling are contiguous, ordered, and start after the left graph's maximum.
- `build_tree_using_neighbor_joining` for 1 and 2 input graphs; clean error for 0.
- Duplicate path name detection: within one graph, across two graphs, and the unnamed-path case.

**Integration** (`packages/pangraph/tests/itest_merge.rs`)

- Split `data/ges-1.fa`, `build` A and B separately, `merge` them, then assert: every original
  sequence is reconstructed byte-identically (matched **by name**, not by order); `sanity_check`
  passes; the path count is `nA + nB`; all path names are preserved.
- A single-genome right graph (exercises phase 1).
- Two graphs with no shared homology — the merge degenerates to a pure join, and this is also the
  case that produces surviving small integer ids, so it is the sharpest test of §3.
- Merging a graph with itself must fail with a duplicate-name error.

---

## 9. Documentation

- This document, plus `dev/design/merge-implementation-notes.md` for as-built decisions.
- A user-facing tutorial section on extending an existing graph, stating the cost model (§5), the
  base-graph-first convention (§3.4), and that `merge(build(A), build(B)) != build(A ∪ B)` (§1).
- `pangraph reconstruct` docs: record order and `index` are only meaningful for graphs produced
  directly by `build` (§4.4).
- Regenerate `docs/docs/reference.md` via `docs/generate-reference-docs` after the CLI changes.
- CHANGELOG entry.

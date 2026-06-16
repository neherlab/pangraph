# Annotation lifting — implementation notes

> Living, **as-built** log for the annotation-lifting feature. The design rationale lives in
> [`annotate.md`](./annotate.md) (the manifesto); this file records concrete decisions, chosen
> dependencies, behaviours and gotchas as each phase lands, so they are not lost when we write the
> user-facing documentation. Update it at the end of every phase.

## Status

| Phase | Scope | State |
|---|---|---|
| P1 | Internal `Feature` model + GFF reader + seqid↔path matching | ✅ done |
| P1.5 | Real-data smoke tests (klebs graph + NCBI GFF annotations) | ✅ done |
| P2 | Inverse coordinate helper `consensus_coords_from_node` (+ round-trip tests) | ✅ done |
| P3 | Node-level lift + `AnnotationWriter` trait (CSV impl) | ✅ done |
| P5.1 | `pangraph annotate` CLI (node-level output, **exact** seqid matching) | ✅ done |
| P4 | Block-level compaction (coordinate agreement) + block CSV writer | ✅ done |
| P5.2 | Wire block-level output into the `annotate` command | ✅ done |
| P5.3 | Docs page + CLI reference regeneration (both output levels) | ⏳ todo |
| P6 | pypangraph consumer + visualization example | ⏳ todo |

> **Ordering note (2026-06-12):** P4 (block compaction) is **deferred until after** a working
> node-level CLI (P5.1), so the still-undecided compaction policy can be designed against real
> node-level output. See the roadmap in [`annotate.md`](./annotate.md) §10.

## Module layout

- `packages/pangraph/src/annotation/` — domain logic
  - `feature.rs` — the format-agnostic `Feature` model + coordinate helper
  - `matching.rs` — `match_features_to_paths`
  - `lift.rs` — `LiftedAnnotation` + `lift_feature` / `lift_features` (P3)
  - `writer.rs` — `AnnotationWriter` trait + `CsvAnnotationWriter` (P3)
- `packages/pangraph/src/pangraph/slice.rs` — `consensus_coords_from_node` (P2), next to its forward
  `interval_node_coords` (within-block consensus↔node coordinate transforms)
- `packages/pangraph/src/io/` — format readers (next to `fasta.rs`)
  - `gff.rs` — `GffReader`
- Fixtures: `data/example.gff`

## Dependencies added

- `noodles` gained the **`gff`** feature (root `Cargo.toml`); resolves to **noodles-gff 0.26.0**
  (via noodles `=0.60.0`). This is the older eager-record API: `gff::Reader::new(bufread).records()`
  yields `io::Result<gff::Record>`.

## P1 decisions & behaviours (these should surface in user docs)

### Coordinate conventions are normalized in the readers
Internal `Feature.interval` is always **0-based, half-open `[start, end)`**.

- **GFF** is 1-based, fully-closed → converted with
  `annotation::feature::interval_from_one_based_inclusive` (`(1, 3)` → `[0, 3)`).

### Strand is `Option<Strand>`
- GFF `+`/`-` → `Some(Forward/Reverse)`; GFF `.`/`?` (unstranded) → `None`.

### GFF specifics
- `id` ← `ID` attribute; `name` ← `Name` attribute; `source` `.` → `None`.
- `records()` automatically skips directives/comments and stops at any trailing `##FASTA` section.
- All attributes are carried into `attributes` (multi-value GFF arrays rendered via their `Display`).
- **Blank-line tolerance (real-world fix):** noodles-gff 0.26 parses a blank line as a record and
  errors (`missing field: Source`), but real NCBI GFF trails a blank line. `read_many` now reads the
  whole input and drops whitespace-only lines before parsing (annotation files are small enough to
  buffer). Found via the klebs smoke test, not the synthetic fixtures.

### Matching (`match_features_to_paths`)
- Builds a path-name → `PathId` lookup **once** (O(features + paths)); paths without a name skipped.
- Resolves each feature's effective name through `seqid_map` (the future `--seqid-map`) before lookup.
- **Aggregates *all* unmatched seqids into a single hard error** (via `make_error!`) rather than
  failing on the first — seqid/path-name mismatch is the #1 user failure mode, so the message lists
  every offending seqid and suggests the seqid map.

## P2 decisions & behaviours (the inverse coordinate helper)

`consensus_coords_from_node(node_coords, edits, block_L) -> (usize, usize)` (`slice.rs`) is the
exact reverse of `interval_node_coords`: walking node → consensus it **adds back** deletion lengths
and **subtracts** insertion lengths for edits preceding each endpoint. Substitutions never move
coordinates.

### Conventions for the non-invertible spots
The forward map is not a bijection (insertion bases have no consensus image; deleted consensus bases
have no node image), so two endpoint cases are resolved by convention:

- **Endpoint strictly inside an insertion** → snap to the insertion's consensus anchor `ins.pos`.
  (Internally the per-endpoint helper also returns an `in_insertion` flag; it is **private** in P2
  and will be promoted in P3 when `LiftedAnnotation` gains per-endpoint flags.)
- **Endpoint on a deletion** (node coordinate where several consensus positions collapse) → exclude
  the gap: the **start** snaps to the deletion's **right edge**, the **end** to its **left edge**.
  A feature that lands entirely inside a deletion therefore yields `start > end` — a signal that it
  collapsed inside a gap (P3 will interpret/flag this).

### Caveats & scope
- **Round-trip is exact only for indel-clean endpoints.** `consensus_coords_from_node ∘
  interval_node_coords == identity` holds when no endpoint lands interior to an indel. The existing
  forward fixtures (`test_node_coords`, `test_interval_node_coords`) deliberately have
  deletion-interior endpoints and are therefore *not* exactly invertible — they are reused only to
  pin the convention result.
- **Within-block only.** This is purely a consensus↔node transform inside one block. The genome→node
  hop (path-offset subtraction, strand flip, circular-wrap modular arithmetic via
  `new_position_circular`) is **not** here — it belongs to the per-feature lift (P3).

## P3 decisions & behaviours (the node-level lift)

`lift_feature(feature, path, graph)` performs the genome→node→consensus transform of §4 of the
manifesto and returns one `LiftedAnnotation` per overlapped node (a *segment*), ordered by genome
coordinate. `lift_features(grouped, graph)` runs it over the `match_features_to_paths` output.

### Coordinate frames in the output
- `node_start`/`node_end` are **consensus-oriented** node-local coordinates (already strand-flipped),
  i.e. exactly the input to `consensus_coords_from_node`; `cons_start`/`cons_end` are block-consensus
  coordinates. All half-open, 0-based.
- Both **nodes** and **features** may wrap the circular origin. Nodes are split by
  `node_coverage_pieces` (`[p0, tot_len)` + `[0, p1)`, or a whole-circle node `p0 == p1`); features
  by `feature_pieces` (see the origin-spanning section below).

### Per-endpoint flags use the **consensus-endpoint frame** (decision)
`start_*` / `end_*` flags refer to the row's `cons_start` / `cons_end` endpoints, **not** the genome
5'/3' ends. On a reverse-strand node the feature's genome-start maps to `cons_end`; `strand_on_consensus`
(the source strand flipped on reverse nodes; `None` stays `None`) lets a consumer map back. (A
companion `feature_strand` field — the *un-flipped* GFF genome strand — was added later by the
crossing-identity fix; see the P4 notes.)
- `start_is_terminus`/`end_is_terminus` — the endpoint is a real feature terminus vs a fragment
  boundary from a node/block split. A fully-covered feature has exactly two real termini (its genome
  start and end); interior segment boundaries are non-termini. (GFF source-`partial` carry-through is
  a later refinement; nodes tile a path with no gaps, so within-graph clipping does not arise.)
- `start_in_insertion`/`end_in_insertion` — the endpoint snapped inside an insertion (no consensus
  image). **Policy: keep-and-flag** (decision), per the manifesto lean. A feature wholly inside an
  insertion is kept with `cons_start == cons_end` and both flags set. (The deletion-collapse case
  `start > end` cannot occur for a real, length>0 feature, since deleted positions have no node bases.)

### `feature_id` / `parent_feature_id`
`parent_feature_id` is the source GFF `ID` (shared across a feature's segments; `None` if absent).
`feature_id` is per-row unique: the parent (or a `"{seqid}:{start}-{end}"` fallback when there is no
`ID`), suffixed `.seg{idx}` for multi-segment features. `n_segments` and `frac_covered` (genome-length
fraction of the feature in this segment) are also emitted.

### Writer
`AnnotationWriter` is a trait (`write_node_annotations(&[LiftedAnnotation])`); the default impl
`CsvAnnotationWriter` writes long-format CSV via the `csv` crate over `create_file_or_stdout` (so `-`
= stdout and `.gz/.bz2/.xz/.zst` compress transparently). **Headers are on** (the shared
`CsvStructFileWriter` forces them off, so the writer uses the `csv` builder directly). `attributes`
renders as one JSON-string column (a JSON array of `[key, value]` pairs, order/duplicate-preserving);
`Option`s render as empty cells; ids as plain numbers; `frac_covered` formatted to 4 decimals.

### Validation
`packages/pangraph/tests/itest_annotate_lift.rs` round-trips on `data/test_graph.json` by
**reconstructing** each genome from the graph (`commands::reconstruct::reconstruct_run::reconstruct`,
no new fixtures), lifting synthetic features (single-node, cross-block, and ones overlapping the
origin-wrapping node in both strands) and asserting the reassembled segment bases equal the genome
substring. Unit tests in `lift.rs` pin exact coordinates for each edge case (strand, deletion,
insertion, in-insertion, multi-segment termini, circular wrap, unstranded).

## P5.1 decisions & behaviours (the `annotate` command)

The `annotate` command (`packages/pangraph/src/commands/annotate/`, mirroring `simplify`) is pure
wiring over the P1–P3 library API — no lift-logic changes. `annotate_run` does: `Pangraph::from_path`
→ for each `--gff` file `GffReader::from_path(..).read_many()` →
`match_features_to_paths(.., &BTreeMap::new())` → `lift_features` →
`CsvAnnotationWriter::new(out, b',').write_node_annotations(..)`.

### CLI surface (intentionally minimal)
- **Graph**: positional `input` (`Option<PathBuf>`, stdin if omitted), like `simplify`.
- **GFF(s)**: `--gff <FILE>...`, `required = true` (≥1), `num_args = 1..` — accepts several paths
  after one flag (`--gff a.gff b.gff`, so shell globs `--gff *.gff` work) and/or the flag repeated
  (`--gff a --gff b`); values accumulate into one `Vec<PathBuf>`. Because the graph is the
  positional arg, give it before the flag (or via stdin) so the variadic does not slurp it as an
  extra GFF; a following flag such as `-o` terminates the list. Transparent decompression by extension.
- **Output**: `-o/--output` (default `-` = stdout); CSV only; compression inferred from extension.
- **No `--seqid-map`, no `--output-level`, no `--format`/`--delimiter`** — deferred to later phases
  (§12 / P4 / P5.2). The CSV delimiter is hard-wired to `,`.

### seqid matching = exact (decision)
Matching passes an **empty** `seqid_map`, so `match_features_to_paths` matches `seqid == path.name`
exactly and aggregates **all** unmatched seqids into one hard error (existing P1 behaviour). On real
NCBI data (versioned `NZ_*.1` seqids vs bare `NZ_*` paths) this **errors by design** — the concrete
signal motivating the §12 relaxation decision. (The P1 matcher error text was reworded to drop the
"provide an explicit seqid-to-path mapping" suggestion — a capability the CLI does not expose yet —
and now simply states that seqids must match the path names.)

### Tests
`packages/pangraph/tests/itest_annotate_cli.rs` drives `annotate_run` end-to-end against
`data/test_graph.json` with a GFF written at test time using a **real** path name (exercising the
`GffReader` path that `itest_annotate_lift.rs` skips): asserts the CSV header + that both features
lift + the genome column; a mismatch case (`data/example.gff`, seqids `chr1`/`chr2`) asserts the loud
error; and a `.csv.gz` output is read back through transparent decompression.

## Origin-spanning features (circular paths)

Real NCBI GFF encodes a feature crossing the replicon origin as a **single record with
`end > sequence length`** (e.g. `5278484..5279188` on a 5,278,493-bp chromosome → wraps to
`[0, 695)`), per
<https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/#origin-spanning-features>.
The P5.1 verification on real klebs data found exactly two such features (the chromosome's origin
gene + its CDS); every other feature was already byte-exact (20,935 / 20,937).

**Behaviour (`lift.rs`):**
- A feature with `f_e > tot_len` on a **circular** path is decomposed by `feature_pieces` into its
  arc pieces `[f_s, tot_len)` + `[0, f_e − tot_len)` (each tagged with an `arc_off`), intersected with
  the node pieces, and emitted as a **single feature** (shared `parent_feature_id`).
- `merge_wrapped_segments` folds consecutive same-node, node-contiguous raw segments back together, so
  a wrap landing on **one** origin-wrapping (or whole-circle) node yields **one** segment — not an
  artificial head/tail split. A genuine multi-block wrap stays multiple segments, ordered 5'→3'.
- Terminus flags and `frac_covered` are computed from **arc position** along the feature
  (`arc_start == 0` / `arc_end == len`), identical to the old genome-coordinate test for non-wrapping
  features.
- **Skipped with a `warn!`** (returns no segments, so the user is told): a wrap on a **non-circular**
  path, a feature **longer than the genome**, or a **start beyond the genome**. Warnings go through the
  `log` crate, surfaced by `--verbosity`.

Re-running the (throwaway) klebs verification after this change: **20,937 / 20,937 OK, 0 skipped**.

## P4 decisions & behaviours (block-level compaction)

`annotation/compact.rs` adds the **block-level** view: collapsing the redundant node-level table into
one **consensus feature** per recurring placement. It is a configurable, **modular** layer built on
top of the node-level lift — never baked into it.

### Modular strategy trait
`BlockCompactionStrategy::compact(&[LiftedAnnotation], &Pangraph) -> Vec<BlockAnnotation>` is the
single entry point. The one P4 impl is `CoordinateConsensusStrategy { min_frequency,
property_threshold }`; future name-/ortholog-based strategies implement the same trait. Provisional
config defaults live on `Default` (`min_frequency = 0.9`, `property_threshold = 0.5`) — the
user-facing defaults are owned by the CLI in **P5.2**.

### Clustering = coordinate-consensus (decisions locked with the user)
> **Superseded by the crossing-identity fix (2026-06-16)** — the two-endpoint key and its
> mixed-orientation consequence below describe the *original* P4 design; the current behaviour is the
> whole-crossing multiset described in **"Crossing-identity fix"** further down. The bullets are kept
> for history.

- **Cluster key** = the feature's two block-consensus **terminus** endpoints **+ `feature_type` +
  `strand_on_consensus`**. A gene vs CDS, or opposite strands, at identical coordinates stay distinct
  (the manifesto §8's coordinate-only key was extended with type + strand on the user's call).
- Each genome's feature instance is recovered by grouping node-level rows on
  `(genome, feature_type, base id)` (`base id` = `parent_feature_id`, else the `feature_id` with its
  `.seg{idx}` suffix stripped), then **reduced to exactly two terminus endpoints** via the
  `*_is_terminus` flags. Instances not yielding exactly two termini (e.g. partial features) are
  **excluded** from compaction — they remain in the node-level table. The two endpoints are stored in
  **canonical sorted order** (`start` = lower `(block, coord)`, not necessarily the genome 5');
  `strand_on_consensus` disambiguates orientation.
- **Consequence:** a block traversed in **mixed orientations** across genomes splits the same
  biological gene into a `+` and a `-` cluster (strand is part of the identity). Acceptable: the
  consensus frame genuinely differs, and most blocks are traversed one way.

### M-of-N threshold (frequency)
- `M` = number of distinct genomes sharing the exact key (a genome counts once even if duplicated).
- `N` = number of distinct genomes **traversing the cluster's block(s)** — the intersection of the
  per-block path sets (just the one block for the common single-block case), via
  `PangraphBlock::isolates`. *Not* the total number of genomes in the graph: a gene is not penalised
  for being absent in genomes that lack the block entirely.
- A cluster is emitted when `M >= max(1, ceil(min_frequency * N))`.

### Consensus metadata
Over the `M` supporters, `consensus_name` is the majority `name` and `consensus_attributes` are the
per-key majority attribute values, each emitted only when its support clears
`ceil(property_threshold * M)` (default 0.5). Ties break deterministically to the smaller string;
output is key-sorted. Resolved generically over **all** attribute keys, so `product` etc. fall out
for free.

### Writer
The `AnnotationWriter` trait gained `write_block_annotations(&[BlockAnnotation])`, implemented on
`CsvAnnotationWriter` with a `BlockAnnotationCsvRow` mirroring the node-level row (`consensus_attributes`
as one JSON-string column, strand as `+`/`-`/empty, ids as numbers, headers on). **CSV only** this
phase; JSON / GFF-on-consensus writers are deferred. CLI wiring is **P5.2** (below).

### Validation
Unit tests in `compact.rs` pin each behaviour (all-agree, below-threshold, type split, strand split,
name/attribute thresholds, multi-block endpoints, `N` = paths-traversing-block, partial-feature
exclusion). `tests/itest_annotate_compact.rs` lifts whole-core-block-node features across every
genome of `data/test_graph.json` and asserts they collapse to consensus feature(s) at `[0, L)` whose
supporters sum to all genomes, then round-trips the block CSV writer.

### Crossing-identity fix (2026-06-16, supersedes the two-endpoint key)
Running `annotate blocks` on real *E. coli* graphs surfaced a bug on the gene `yjeM` (see the
investigation in the PR): a feature whose body sits on a core block and whose 37 bp tail sits on an
**inverted** block was split into separate `+` and `-` clusters even though the placement was
*identical* across genomes — and one genome was dropped below threshold. Root cause: the cluster key
sampled a single representative strand from `segment_idx 0` (the genome-lowest segment), and that
segment flips with the genome's global orientation across an inversion boundary. Strand was being
treated as one per-feature value when it is really **per-segment**.

What changed:
- **Identity = the whole crossing.** The cluster key is now `(feature_type, Vec<Segment>)` where
  `Segment = (BlockId, cons_start, cons_end, strand_on_consensus)`, ordered **5'→3'**. Per-block
  `strand_on_consensus` is invariant across genomes (the genome's `+/-` is absorbed by node
  orientation), so the whole-crossing multiset is stable. This also distinguishes `X→Y→X` from `X→Y`.
- **Ordering** uses the new node-level `feature_strand` field (the original GFF genome strand,
  un-flipped) to reverse the genome-ordered segments for reverse-strand features; unstranded features
  (no reading direction) are canonicalized by the lexicographically smaller of the two orders. We
  chose to **add `feature_strand`** rather than recompute `strand_on_consensus XOR node.strand()` at
  compaction, because it keeps compaction independent of graph-consistent `node_id`s (its unit tests
  build `LiftedAnnotation`s directly) and makes the node-level CSV self-document the GFF strand.
- **Output is one row per segment**, sharing a `cluster_id` and numbered 5'→3' by `segment_idx`
  (`n_segments` total). A **duplicated block crossed twice** keeps both occurrences as separate rows,
  distinguished by `segment_idx`.
- **`N` is structural capability** — genomes traversing every block in the crossing with the required
  **multiplicity**, where multiplicity is the crossing's count of **distinct nodes** per block (not
  segments), reduced to the element-wise minimum over the supporters so `N >= M` always holds.
  Counting by node, not segment, matters: a whole-genome `region` feature crosses an **origin-spanning
  node** that the lift splits into two coverage pieces — two segments on **one** node — so a
  segment-count requirement of 2 against a genome with one such node gave `N = 0` while `M = 1` (an
  `M > N` contract violation found on the E. coli data). `block_genome_counts` supplies per-genome
  node counts; `reduce_instance` records per-block distinct-node use. A body-shared/tail-different
  feature thus fragments into several clusters, each scored only against the genomes that could carry
  it, so every variant is promoted at its true support.
- **Per-segment support** `n_support_segment / n_total_segment` is reported per row, computed over the
  whole node-level table independently of clustering (`segment_support_sets`), so a body segment
  shared by several crossings reports the same value in each. It is informational; gating stays
  per-crossing.
- New unit tests: `test_inversion_crossing_collapses_to_one_cluster` (the regression),
  `test_duplicated_block_crossed_twice`, `test_shared_body_three_tails_each_confident`; the existing
  tests were updated to the per-segment shape and the block CSV header is now
  `type,cluster_id,segment_idx,…`.

### Feature-type filtering (`--only-type` / `--exclude-type`)
Motivated by the same E. coli run: whole-contig `region` records (one per genome) expanded to ~30k of
~65k block rows — correct output, but noise. Rather than special-case any type in the compactor, the
input is filtered by GFF `type` up front. Two **shared** `AnnotateCommonArgs` flags
(`packages/pangraph/src/commands/annotate/annotate_args.rs`), so both `nodes` and `blocks` inherit
them: `--only-type` (whitelist) and `--exclude-type` (blacklist), each `value_delimiter = ','`
(comma-separated and/or repeatable, accumulating into a `Vec<String>`) and `conflicts_with` each other
(mutually exclusive). Matching is **exact and case-sensitive**. The filter — `filter_features_by_type`
in `annotation/feature.rs`, a pure `Vec<Feature> -> Vec<Feature>` — runs in `load_and_lift` **after**
GFF reading but **before** `match_features_to_paths`, so an excluded type never triggers a
seqid-mismatch error and never reaches the lift. Empty lists (the default) are a no-op, so existing
behaviour is unchanged. Unknown type names are silently no-ops (whitelisting an absent type yields an
empty result; excluding one keeps everything).

## P5.2 decisions & behaviours (block output on the CLI)

P5.2 exposes compaction on the `annotate` command. Following the `export` precedent (one command,
several output shapes as subcommands), `annotate` became a **subcommand group**:

- `pangraph annotate nodes` — the P5.1 node-level table (unchanged behaviour).
- `pangraph annotate blocks` — runs `CoordinateConsensusStrategy::compact` then
  `write_block_annotations`.

### CLI surface (decisions locked with the user)
- **Subcommands, not a `--output-level` flag** — keeps block-only flags out of node-mode `--help`,
  and matches `export`. Bare `annotate` requires a subcommand (no default mode; acceptable while
  `annotate` is unreleased on the integration branch).
- **Names** = `nodes` / `blocks` (plural, consistent).
- **Shared options** (`input`, `--gff`, `-o/--output`) live in a flattened `AnnotateCommonArgs`;
  both subcommands embed it via `#[clap(flatten)]` (DRY vs `export`'s per-variant duplication).
- **Block-only flags**: `--min-frequency` (default 0.9) and `--property-threshold` (default 0.5),
  whose defaults are kept in sync with `CoordinateConsensusStrategy::default()`. Both are validated at
  parse time by a `parse_fraction` `value_parser` rejecting anything outside `[0, 1]` (incl. `NaN`/
  `inf`) with a clap error, rather than silently emitting nothing / promoting everything downstream.
  A `--strategy` selector is **deferred** until a second strategy exists — a one-option flag is noise,
  and the trait already provides the modularity internally.

### Wiring
`annotate_run` matches the `PangraphAnnotateArgs` enum and dispatches to `annotate_run_nodes` /
`annotate_run_blocks`. Both share `load_and_lift(&AnnotateCommonArgs) -> (Pangraph, Vec<LiftedAnnotation>)`
(graph load → per-`--gff` `read_many` → `match_features_to_paths` → `lift_features`); the graph is
returned because block compaction needs it for per-cluster path totals. `tests/itest_annotate_cli.rs`
exercises both subcommands (`annotate blocks --min-frequency 0` emits the block-level header + rows)
and pins the threshold wiring: with a single annotated genome, `--min-frequency 1.0` drops the
single-support clusters that `0.0` keeps (a field swap would surface as the lenient run coming back
empty).

## Public API introduced in P1–P5.2

```rust
// annotation::feature
pub struct Feature { seqid, source, feature_type, interval, strand: Option<Strand>, id, name, attributes }
pub fn interval_from_one_based_inclusive(start: usize, end: usize) -> Interval;

// io::gff
pub struct GffReader<'a>;          // new / from_str / from_path / read_many() -> Vec<Feature>

// annotation::matching
pub fn match_features_to_paths(
  features: Vec<Feature>, graph: &Pangraph, seqid_map: &BTreeMap<String, String>,
) -> Result<BTreeMap<PathId, Vec<Feature>>, Report>;

// pangraph::slice (P2) — inverse of interval_node_coords
pub fn consensus_coords_from_node(node_coords: (usize, usize), edits: &Edit, block_L: usize) -> (usize, usize);

// pangraph::slice (P3) — same map, exposing per-endpoint in-insertion flags
pub fn consensus_coords_from_node_flagged(
  node_coords: (usize, usize), edits: &Edit, block_L: usize,
) -> ((usize, bool), (usize, bool));

// annotation::lift (P3)
pub struct LiftedAnnotation { /* feature_id, parent_feature_id, segment_idx, n_segments, genome,
  block_id, node_id, strand_on_consensus, node_start/_end, cons_start/_end, start/end_is_terminus,
  start/end_in_insertion, frac_covered, feature_type, name, attributes */ }
pub fn lift_feature(feature: &Feature, path: &PangraphPath, graph: &Pangraph) -> Result<Vec<LiftedAnnotation>, Report>;
pub fn lift_features(grouped: &BTreeMap<PathId, Vec<Feature>>, graph: &Pangraph) -> Result<Vec<LiftedAnnotation>, Report>;

// annotation::compact (P4)
pub struct BlockAnnotation { /* feature_type, strand_on_consensus, start_block_id, cons_start,
  end_block_id, cons_end, consensus_name, consensus_attributes, n_support (M), n_total (N) */ }
pub trait BlockCompactionStrategy {
  fn compact(&self, node_annotations: &[LiftedAnnotation], graph: &Pangraph) -> Result<Vec<BlockAnnotation>, Report>;
}
pub struct CoordinateConsensusStrategy { pub min_frequency: f64, pub property_threshold: f64 } // + Default

// annotation::writer (P3 trait, P4 block method)
pub trait AnnotationWriter {
  fn write_node_annotations(&mut self, annotations: &[LiftedAnnotation]) -> Result<(), Report>;
  fn write_block_annotations(&mut self, annotations: &[BlockAnnotation]) -> Result<(), Report>; // P4
}
pub struct CsvAnnotationWriter;    // new(filepath, delimiter) ; the default CSV impl (node + block)

// commands::annotate::annotate_args (P5.1 → P5.2)
pub enum PangraphAnnotateArgs { Nodes(PangraphAnnotateNodesArgs), Blocks(PangraphAnnotateBlocksArgs) }
pub struct AnnotateCommonArgs { input: Option<PathBuf>, gff: Vec<PathBuf>, output: PathBuf }
pub struct PangraphAnnotateNodesArgs { common: AnnotateCommonArgs }
pub struct PangraphAnnotateBlocksArgs { common: AnnotateCommonArgs, min_frequency: f64, property_threshold: f64 }
```

## Real-data smoke tests & fixtures (P1.5)

`packages/pangraph/tests/itest_klebs_annotations.rs` exercises the readers + matcher on real NCBI
data:

- **`data/klebs_graph.json.gz`** — pangraph built from all 9 genomes in `data/klebs.fa.gz`
  (`pangraph build -c`, ~2m20s, 9 paths / 1378 blocks). Paths are named by the **bare** RefSeq
  accession (the FASTA record id), e.g. `NZ_CP013711`.
- **`data/klebs_annotations/{NZ_CP013711,NC_017540}.gff.gz`** — RefSeq GFF annotations for two of
  those genomes (NCBI sviewer `report=gff3`). Only 2 kept to bound fixture size.
- Tests: parse each GFF (thousands of features, CDS/strand/name present) and **match + lift both
  against the graph** by **exact** seqid equality (no seqid map), asserting the lifted segments stay
  in-bounds (P5.1 real-data lift smoke).

**Key real-world finding — seqid version mismatch.** As downloaded, annotation seqids are the
**versioned** accession (`NZ_CP013711.1`) while graph path names are the **bare** accession
(`NZ_CP013711`). The committed fixtures have since been **normalized** (the `.N` version stripped from
the seqid column only, leaving attribute values intact), so the smoke test matches by exact equality
without a map. The underlying NCBI reality still holds for *user* data, so it remains the motivation
for the §12 version-insensitive-matching decision — it is the default situation for NCBI downloads.

## Carried-forward items for later phases / user docs

> Several of these are now consolidated into the **post-prototype, pre-docs punch list** in
> [`annotate.md`](./annotate.md) §12 (notably seqid↔path matching and the sequence-version-drift
> guard) — settle them there before P5.3.

- Document the `--seqid-map` file format (P5) and the "seqids must match FASTA record names" rule.
- **Version-insensitive seqid matching (P5):** strongly consider auto-stripping the `.N` version so
  `NZ_CP013711.1` matches a `NZ_CP013711` path without an explicit map (see finding above).
- **Sequence-version drift affects the lift (P2/P3), not matching:** when an annotation's source
  record differs in length from the graph's genome (observed with the earlier russian-doll plasmid
  attempt: e.g. `NZ_CP011582` 43433 bp vs graph 45279 bp), coordinates past the divergence lift
  incorrectly. Annotations should come from the *same* sequence used to build the graph; worth a
  user-doc warning and possibly a length/identity sanity check.
- Decide handling of features wholly inside an insertion relative to consensus (P3 open question).
- Note unstranded-GFF handling in user-facing docs.

## GenBank support deferred (v1 decision)

v1 ships **GFF only**; the GenBank reader (the `gb-io` dependency, ~225 lines, and the `.gbk`
fixtures) was **removed** after an initial implementation. Why:

- The lift uses the **graph's** sequence, so GenBank's bundled sequence is unused in this workflow.
- GenBank's compound-location model (`join`/`order`, fuzzy ends) is lossy precisely where the lift is
  sensitive: collapsing to the outer span mis-places an origin-spanning `join(…,1..N)` on a circular
  replicon (≈ the whole genome). GFF encodes multi-segment features as separate clean records — a
  better fit for the per-segment `Feature` model.
- The `.gbk.gz` fixtures dominated the annotation fixture footprint (~7 MB).

**Re-adding it later is cheap** thanks to the format-neutral `Feature`: add a new `io/genbank.rs`
reader producing `Vec<Feature>`, re-add `gb-io`, and handle compound/origin-wrap locations properly
(not a naive outer-span collapse). **User workaround:** convert GenBank → GFF3 with `bp_genbank2gff3`,
EMBOSS `seqret`, or by re-exporting from bakta/prokka.

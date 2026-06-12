# Design manifesto: lifting genome annotations onto the pangenome graph

> Status: **draft / design only** — no implementation yet. This document is the reference we
> execute against incrementally. It is intentionally internal (not published to docs.pangraph.org).

## 1. Motivation & use cases

Most genomes fed into PanGraph come with annotations (genes, CDS, features), most commonly as
**GFF** files. Analysts want to study those annotations *on the graph* rather than on each
isolated genome:

- Which genes fall in **diverse / accessory** regions vs the conserved core?
- Where are the **resistance genes**, and how is their neighbourhood rearranged across genomes?
- **Visualizing** features on a block's consensus (an annotation track per block).

The primary use cases are **visualization** and **locating interesting features on the genomes /
graph**. There is currently no way to express an annotation — defined in genome coordinates — in
**graph coordinates**.

Scope of the feature: a new `pangraph annotate` command that takes an existing graph plus a set of
**GFF** files and produces **lifted annotations**, in two flavours:

- **Node-level** — the lossless source of truth: every feature instance placed on the node(s) it
  overlaps, with consensus coordinates.
- **Block-level (compacted)** — an aggregated, consensus-centric view: per block, the features
  present, with a summary of where they sit on the consensus and **whether the underlying genomes
  agree** on those coordinates.

Two sub-problems, of increasing difficulty:

1. **Map annotations onto graph nodes** — *easy*. A node stores its explicit interval on the path,
   so intersecting a feature interval with node intervals is direct.
2. **Assign annotations to a block's consensus** — *harder*. A node's sequence differs from the
   block consensus by edits (substitutions, insertions, deletions) and possibly strand, so feature
   coordinates must be translated through the alignment.

## 2. Overarching design principle: stay implementation-agnostic

Precise formats and implementations will change as we iterate. Decouple three concerns so any one
can be swapped without touching the others:

- **Parsing (input)** — normalize the input (GFF) into **one internal `Feature` model**. The parser
  is an implementation detail behind that model, and the model is deliberately format-neutral so
  other readers (e.g. GenBank) can be added later without touching the lift (see §11).
- **The lifted result (core)** — first-class in-memory types, e.g. `LiftedAnnotation` (node-level)
  and `BlockAnnotation` (block-level), carrying coordinates + flags + flexible metadata. The lift
  logic produces *these objects*, never format-specific rows.
- **Serialization (output)** — a pluggable **`AnnotationWriter` trait** with concrete impls (CSV,
  JSON, and later e.g. GFF-on-consensus or BED). Adding a format = adding a writer; the lift does
  not change. **CSV is the default** impl; **JSON is one alternate** impl.

```text
GFF ──parse──▶ Feature (internal)
                              │
                              ▼  lift (coordinate transform)
                  LiftedAnnotation / BlockAnnotation
                              │
                              ▼  AnnotationWriter (trait)
                   CSV (default) │ JSON (opt) │ future formats
```

## 3. Coordinate model

Three coordinate frames, traversed in order:

1. **Genome / path coordinate** — what the GFF feature is defined in.
2. **Node-local coordinate** — position within a single node's slice of that genome. Genome→node is
   a direct subtraction of the node's path offset (`node.position.0`), plus a strand flip and modular
   arithmetic for circular paths.
3. **Block-consensus coordinate** — node-local→consensus, obtained by inverting the node's edits.

The natural, compact **graph coordinate** of a feature is therefore:

```
(block_id, cons_start, cons_end, strand_on_consensus)
```

This is an advantage of PanGraph's block model: each block's consensus is a ready-made per-block
reference, so a lifted feature has a clean, compact address — unlike base-level graphs.

## 4. The lift algorithm (per feature)

Given a feature on path `P` with interval `[f_s, f_e)` (0-based half-open) and strand:

1. **Find overlapping nodes.** Iterate the nodes of path `P`; each `PangraphNode` stores
   `position: (usize, usize)` on the path. Intersect `[f_s, f_e)` with each node interval. A feature
   overlapping several nodes **splits into segments** sharing a `parent_feature_id` with a
   `segment_idx`. Handle circular wrap (a node or feature may cross the origin).
2. **Genome → node-local.** Subtract `node.position.0` (modular for circular paths) to get
   `[n_s, n_e)` within the node's genome-oriented sequence.
3. **Strand.** If `node.strand` is reverse, the genome-oriented sequence is the reverse-complement of
   the consensus-oriented node sequence `S_node` (length `len_node = position.1 - position.0`), so
   convert offsets `n_s' = len_node - n_e`, `n_e' = len_node - n_s` and flip the feature's strand.
   Forward strand: unchanged.
4. **Node-seq → consensus.** Invert the edits: as we go from node coordinates back to consensus
   coordinates, **add back deletion lengths** and **subtract insertion lengths** for edits preceding
   the position. This is the **inverse of `interval_node_coords`** (see §6). An endpoint that lands
   strictly inside an insertion has **no consensus position** → snap it to the insertion anchor
   (`ins.pos`) and set an `in_insertion` flag.
5. **Emit** a `LiftedAnnotation` per segment, with consensus coordinates and flags.

### Per-endpoint terminus flags

Each node-level annotation records, **per endpoint independently**, whether the coordinate is the
feature's *real* terminus or merely a *fragment* boundary introduced by a block/node split:
`start_is_terminus` and `end_is_terminus`. This is distinct from `segment_idx`. Example: a gene
split across 3 blocks yields 3 segments; only the **first segment's start** and the **last segment's
end** are true termini — every interior boundary is a fragment end. (Genome features already flagged
partial in the source GFF are carried through as non-termini too.)

## 5. Edge cases (each gets an explicit flag/field)

- **Insertions** — node bases absent from the consensus have no consensus coordinate (`in_insertion`).
- **Deletions** — consensus positions absent from the node are skipped by the inverse map.
- **Substitutions** — do not move coordinates (identity differs, position does not); irrelevant to the lift.
- **Reverse strand** — offset flip + feature-strand flip (§4 step 3).
- **Circular wrap** — feature and/or node crossing the origin; use `PangraphPath.tot_len` +
  `circular` and the modular pattern from `new_position_circular`.
- **Multi-block features** — split into segments with `parent_feature_id` + `segment_idx`; interior
  boundaries are fragment ends.
- **Partial / truncated features** — recorded via the per-endpoint terminus flags and `frac_covered`.

## 6. Reuse map — existing code to build on (do not reinvent)

The hard part is largely solved by existing infrastructure:

- **`interval_node_coords(i, edits, block_L)`** — `packages/pangraph/src/pangraph/slice.rs:103`.
  Maps a **consensus interval → node-sequence coordinates** by walking deletion/insertion lengths.
  Our lift is its **inverse** (node-seq → consensus). Mirror this function to add
  `consensus_coords_from_node(interval, edits, block_L)` and unit-test the two against each other.
- **`Edit` / `Sub` / `Del` / `Ins`** — `packages/pangraph/src/pangraph/edits.rs:17-120`. All edit
  positions are in **consensus coordinates**. Each provides `reverse_complement(len)` — reuse for
  strand handling instead of re-deriving.
- **`Edit::apply`** — `packages/pangraph/src/pangraph/edits.rs:307`. Forward transform
  (consensus → node sequence); used in round-trip validation (§7).
- **`PangraphNode`** — `packages/pangraph/src/pangraph/pangraph_node.rs:19-25`. Explicit
  `position: (usize, usize)` + `strand`; `start_is_end()` (`:59`) flags the circular whole-path node.
- **`new_position_circular` / `new_position_non_circular`** — `slice.rs:67-101`. Strand- and
  circularity-aware position mapping patterns to mirror.
- **`PangraphPath`** — `tot_len` + `circular` for wrap arithmetic; `name` for seqid matching.

### Command wiring & IO reuse

- **Command pattern** (mirror the simplest, `simplify`, in `commands/simplify/`): `{cmd}_args.rs`
  (clap `Parser` struct) + `{cmd}_run.rs` (`pub fn {cmd}_run(args) -> Result<(), Report>`), then
  register the variant in `commands/root_args.rs` (`PangraphCommands` enum) and dispatch in
  `commands/main.rs`.
- **CSV writer** — `CsvStructFileWriter` (`packages/pangraph/src/io/csv.rs:31`) writes serde rows;
  back the default `AnnotationWriter` impl with it.
- **JSON / graph load** — `json_write_file` (`io/json.rs`), `Pangraph::from_path`.
- **Transparent compression** — `create_file_or_stdout` / `open_file_or_stdin`
  (`packages/pangraph/src/io/file.rs:45-80`) auto-handle gz/bz2/xz/zst by extension, for both
  reading and writing, and support `-` for stdin/stdout. Output compression is free.
- **Parsing crates** — `noodles` is already a dependency (its `sam` feature is used in
  `align/bam/cigar.rs`); enable its **`gff`** feature for GFF3. Normalize into the internal `Feature`
  model immediately. (GenBank parsing is deferred for v1 — see §11.)
- **Tests** — integration tests live in `packages/pangraph/tests/itest_*.rs`, using `rstest`
  `#[case]`, `tempfile::tempdir()`, and `pretty_assertions`, against `data/test_graph.json`.

## 7. Command interface

`pangraph annotate`:

- **Inputs**: a graph JSON + one or more GFF files.
- **Mode**: select node-level vs block-level output.
- **Output**: path (default `-` = stdout); compression inferred from the file extension; format
  chosen by the `AnnotationWriter` (CSV default, JSON optional).
- **`--seqid-map`** escape hatch: map a GFF `seqid` to a PanGraph path name.
  **ID matching is the #1 failure mode** — fail loudly (error, not silent drop) on any annotation
  whose seqid does not match a path.

Wire it by mirroring `simplify` (args/run), then register in `root_args.rs` + `main.rs`.

## 8. Annotation types + output schemas

Writers serialize the in-memory types; **schema ≠ format**.

### Node-level annotation (`LiftedAnnotation`) — default, lossless source of truth

Fields:

```
feature_id, parent_feature_id, segment_idx,
genome (path name), block_id, node_id,
strand_on_consensus,
node_start, node_end,            # node-local coordinates
cons_start, cons_end,            # block-consensus coordinates
start_is_terminus, end_is_terminus,   # real terminus vs fragment boundary (per endpoint)
in_insertion,                    # endpoint(s) fell inside an insertion
frac_covered,                    # fraction of the original feature retained in this segment
type, name,                      # e.g. CDS / gene, and gene name
attributes                       # flexible key→value metadata
```

`attributes` keeps us flexible on metadata: the **CSV writer** serializes it as a JSON string in
one column; the **JSON writer** keeps it nested. Same `LiftedAnnotation` objects, different
`AnnotationWriter` impl — the CSV writer renders one row per annotation, the JSON writer renders
them nested.

### Block-level annotation (`BlockAnnotation`) — compacted

One entry per (block, clustered feature). Carries a consensus-coordinate summary plus the
**coordinate-agreement** information across the genomes that carry the feature:

```
block_id, feature (cluster key: name/product),
cons_start_median, cons_start_min, cons_start_max,
cons_end_median,   cons_end_min,   cons_end_max,
n_support,         # genomes carrying this feature on this block
n_total,           # genomes traversing this block
coords_agree       # all supporting genomes agree on start/end (easy case) vs disagree
```

Default output is **long-format CSV**; an **optional nested JSON** carries the full per-genome
coordinate distribution. Clustering "the same annotation" across genomes (by name/product +
consensus overlap) is **opinionated** — keep it a clearly documented, configurable layer **on top
of** the node-level lift, never baked into it.

## 9. Validation strategy

Mirror the build pipeline's "verify by reconstruction" philosophy:

- **Round-trip per feature**: lift a feature to consensus coordinates, push it back through
  `Edit::apply` + strand + node offset, extract that genome subsequence, and compare to the gene
  sequence read from the input FASTA used to build the graph. A match proves the lift correct — this catches every
  strand / indel / wrap bug cheaply.
- **Unit tests** (`rstest`) for the inverse coordinate helper `consensus_coords_from_node` on
  hand-built edits: substitution-only, deletion, insertion, insertion-at-boundary, reverse strand,
  and circular-wrap cases — asserting it inverts `interval_node_coords`.

## 10. Phased implementation roadmap

Land incrementally, one PR/commit per phase:

- **P1** — internal `Feature` model + GFF reader (noodles `gff`) + seqid↔path matching (with
  `--seqid-map`, hard error on mismatch). *(GenBank reader deferred — see §11.)*
- **P2** ✅ — inverse coordinate helper `consensus_coords_from_node` in `slice.rs` + unit tests
  (round-trip against `interval_node_coords`).
- **P3** — per-feature node-level lift (overlap, segmentation, strand, per-endpoint terminus flags)
  producing `LiftedAnnotation` objects + the `AnnotationWriter` trait with a CSV impl + integration
  test on `data/test_graph.json`.
- **P4** — block-level compaction (clustering + coordinate-agreement) into `BlockAnnotation` objects
  + a second writer impl (CSV default, JSON optional) over the same trait + tests.
- **P5** — `annotate` command wiring (args/run/register/dispatch), Docusaurus docs page, CLI
  reference regeneration.
- **P6 (future)** — pypangraph consumer (load the tables; feature→node→block→path joins) + a
  visualization example. Out of scope for this Rust manifesto beyond this forward pointer.

## 11. Open questions / decisions deferred

- **Input formats (v1 = GFF only)** — v1 supports **GFF** only; a GenBank reader is **deferred**.
  Rationale: the lift uses the *graph's* sequence, so GenBank's one real advantage — bundling
  sequence + annotation — is unused here; GenBank's compound-location model (`join`/`order`, fuzzy
  ends) is lossy exactly where the lift is sensitive (an origin-spanning `join(…,1..N)` on a circular
  replicon collapses to ≈ the whole genome), whereas GFF encodes multi-segment features as separate
  clean records that fit the per-segment `Feature` model better; and the `.gbk` fixtures are large.
  Because `Feature` is format-neutral, re-adding a GenBank reader later is a self-contained job (a new
  `io/genbank.rs` producing `Vec<Feature>`, handling compound/origin-wrap properly rather than
  collapsing to the outer span). Users who only have GenBank can convert to GFF3 (`bp_genbank2gff3`,
  EMBOSS `seqret`, or re-export from bakta/prokka).
- **Embed in graph JSON vs separate file** — default to a separate file; revisit if a single
  self-contained artifact is wanted.
- **Block-level clustering policy** — by feature name/product, by consensus overlap, or by an
  external ortholog grouping; make it configurable.
- **Features wholly inside an insertion** — drop, or keep node-level-only with no consensus
  coordinate? Lean towards keep-and-flag.
- **Coordinate convention** — 0-based half-open internally; convert only at GFF I/O (GFF is 1-based
  closed) to avoid off-by-one creep.

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
| P4 | Block-level compaction (coordinate agreement) + JSON writer | ⏳ todo |
| P5 | `pangraph annotate` CLI command + `--seqid-map` + docs/CLI reference | ⏳ todo |
| P6 | pypangraph consumer + visualization example | ⏳ todo |

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
- A feature interval is assumed **non-wrapping** (`f_s < f_e`); origin-spanning features are separate
  GFF records. **Nodes**, however, may wrap the circular origin and that is handled
  (`node_coverage_pieces` splits a wrapping node into its `[p0, tot_len)` and `[0, p1)` pieces, and a
  whole-circle node `p0 == p1`).

### Per-endpoint flags use the **consensus-endpoint frame** (decision)
`start_*` / `end_*` flags refer to the row's `cons_start` / `cons_end` endpoints, **not** the genome
5'/3' ends. On a reverse-strand node the feature's genome-start maps to `cons_end`; `strand_on_consensus`
(the source strand flipped on reverse nodes; `None` stays `None`) lets a consumer map back.
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

## Public API introduced in P1–P3

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

// annotation::writer (P3)
pub trait AnnotationWriter { fn write_node_annotations(&mut self, annotations: &[LiftedAnnotation]) -> Result<(), Report>; }
pub struct CsvAnnotationWriter;    // new(filepath, delimiter) ; the default CSV impl
```

## Real-data smoke tests & fixtures (P1.5)

`packages/pangraph/tests/itest_klebs_annotations.rs` exercises the readers + matcher on real NCBI
data:

- **`data/klebs_graph.json.gz`** — pangraph built from all 9 genomes in `data/klebs.fa.gz`
  (`pangraph build -c`, ~2m20s, 9 paths / 1378 blocks). Paths are named by the **bare** RefSeq
  accession (the FASTA record id), e.g. `NZ_CP013711`.
- **`data/klebs_annotations/{NZ_CP013711,NC_017540}.gff.gz`** — RefSeq GFF annotations for two of
  those genomes (NCBI sviewer `report=gff3`). Only 2 kept to bound fixture size.
- Tests: parse each GFF (thousands of features, CDS/strand/name present) and **match both against
  the graph**, building the seqid map from the files themselves.

**Key real-world finding — seqid version mismatch.** Annotation seqids are the **versioned**
accession (`NZ_CP013711.1`); graph path names are the **bare** accession (`NZ_CP013711`). The smoke
test bridges this with a `version → bare` seqid map. This strongly motivates **version-insensitive
matching** (or a documented `--seqid-map`) in P5, since it is the default situation for NCBI data.

## Carried-forward items for later phases / user docs
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

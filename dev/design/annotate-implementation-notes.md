# Annotation lifting — implementation notes

> Living, **as-built** log for the annotation-lifting feature. The design rationale lives in
> [`annotate.md`](./annotate.md) (the manifesto); this file records concrete decisions, chosen
> dependencies, behaviours and gotchas as each phase lands, so they are not lost when we write the
> user-facing documentation. Update it at the end of every phase.

## Status

| Phase | Scope | State |
|---|---|---|
| P1 | Internal `Feature` model + GFF/GenBank readers + seqid↔path matching | ✅ done |
| P1.5 | Real-data smoke tests (klebs graph + NCBI annotations) | ✅ done |
| P2 | Inverse coordinate helper `consensus_coords_from_node` (+ round-trip tests) | ⏳ todo |
| P3 | Node-level lift + `AnnotationWriter` trait (CSV impl) | ⏳ todo |
| P4 | Block-level compaction (coordinate agreement) + JSON writer | ⏳ todo |
| P5 | `pangraph annotate` CLI command + `--seqid-map` + docs/CLI reference | ⏳ todo |
| P6 | pypangraph consumer + visualization example | ⏳ todo |

## Module layout

- `packages/pangraph/src/annotation/` — domain logic
  - `feature.rs` — the format-agnostic `Feature` model + coordinate helper
  - `matching.rs` — `match_features_to_paths`
  - *(later: `lift.rs`, `writer.rs`)*
- `packages/pangraph/src/io/` — format readers (next to `fasta.rs`)
  - `gff.rs` — `GffReader`
  - `genbank.rs` — `GenbankReader`
- Fixtures: `data/example.gff`, `data/example.gbk`

## Dependencies added

- `noodles` gained the **`gff`** feature (root `Cargo.toml`); resolves to **noodles-gff 0.26.0**
  (via noodles `=0.60.0`). This is the older eager-record API: `gff::Reader::new(bufread).records()`
  yields `io::Result<gff::Record>`.
- **`gb-io = "=0.9.0"`** for GenBank: `SeqReader::new(read)` iterates `Result<Seq, _>`.

## P1 decisions & behaviours (these should surface in user docs)

### Coordinate conventions are normalized in the readers
Internal `Feature.interval` is always **0-based, half-open `[start, end)`**.

- **GFF** is 1-based, fully-closed → converted with
  `annotation::feature::interval_from_one_based_inclusive` (`(1, 3)` → `[0, 3)`).
- **gb-io is already 0-based half-open** — a single GenBank position `1` is `Range((0,..),(1,..))`,
  so GenBank needs **no** conversion. (Easy to get wrong — gb-io does the 1-based→0-based shift for us.)

### Strand is `Option<Strand>`
- GFF `+`/`-` → `Some(Forward/Reverse)`; GFF `.`/`?` (unstranded) → `None`.
- GenBank features are always stranded → `Some(...)`; `complement(...)` location → `Reverse`,
  otherwise `Forward`.

### GenBank specifics
- **seqid precedence:** `VERSION` (accession.version) → `ACCESSION` → `LOCUS` name. Mismatches with
  FASTA/path names are expected and handled by the (future) `--seqid-map` override.
- **Compound locations** (`join`/`order`) are collapsed to their **outer span** via
  `Location::find_bounds()` — one `Feature` per GenBank feature. Adequate for bacterial genomes
  (no introns); documented limitation, revisit if spliced features matter.
- Features whose location bounds cannot be resolved are **skipped with a `warn!`** (not a hard error).
- All features are emitted, including the whole-record `source` feature; callers filter by
  `feature_type` if they want only genes/CDS.
- `id` ← `locus_tag`; `name` ← `gene` then falls back to `product`; all qualifiers are kept in
  `attributes` (order- and duplicate-preserving `Vec<(String, String)>`).

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

## Public API introduced in P1

```rust
// annotation::feature
pub struct Feature { seqid, source, feature_type, interval, strand: Option<Strand>, id, name, attributes }
pub fn interval_from_one_based_inclusive(start: usize, end: usize) -> Interval;

// io::gff
pub struct GffReader<'a>;          // new / from_str / from_path / read_many() -> Vec<Feature>

// io::genbank
pub struct GenbankReader<'a>;      // new / from_str / from_path / read_many() -> Vec<Feature>

// annotation::matching
pub fn match_features_to_paths(
  features: Vec<Feature>, graph: &Pangraph, seqid_map: &BTreeMap<String, String>,
) -> Result<BTreeMap<PathId, Vec<Feature>>, Report>;
```

## Real-data smoke tests & fixtures (P1.5)

`packages/pangraph/tests/itest_klebs_annotations.rs` exercises the readers + matcher on real NCBI
data:

- **`data/klebs_graph.json.gz`** — pangraph built from all 9 genomes in `data/klebs.fa.gz`
  (`pangraph build -c`, ~2m20s, 9 paths / 1378 blocks). Paths are named by the **bare** RefSeq
  accession (the FASTA record id), e.g. `NZ_CP013711`.
- **`data/klebs_annotations/{NZ_CP013711,NC_017540}.{gff,gbk}.gz`** — RefSeq annotations for two of
  those genomes (GFF via NCBI sviewer `report=gff3`; GenBank via efetch `gbwithparts`). Only 2 kept
  to bound fixture size (~8 MB).
- Tests: parse each GFF/GenBank (thousands of features, CDS/strand/name present) and **match both
  against the graph**, building the seqid map from the files themselves.

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
- Note the GenBank compound-location collapse and unstranded-GFF handling in user-facing docs.

# Block-compaction policy: open questions to settle before it's final

> Scope: the **block-level compaction** layer only (`annotation/compact.rs`,
> `CoordinateConsensusStrategy`). The **node-level lift is not affected** by any of this — it is the
> lossless source of truth, and compaction is an opinionated summary built on top of it behind the
> `BlockCompactionStrategy` trait. Everything below can change (or a whole alternative strategy can be
> added) without touching the lift or the writers. Resolve these before P5.3 freezes the behaviour in
> user docs. Companion to [`annotate.md`](./annotate.md) §8/§11/§12 and the as-built notes in
> [`annotate-implementation-notes.md`](./annotate-implementation-notes.md).

The layer already works and is regression-tested against the hard structural cases (repeated blocks,
inversions, origin-spanning nodes). What is still *provisional* is the **reconciliation policy**: how
strictly two genomes' placements must agree to be called "the same feature", and what the block-level
output should carry. These are decisions, not bugs.

## 1. Coordinate-exact identity vs. a tolerance (the crux)

> **Resolved (2026-07-06): ship exact-only for v1.** The five-species fragmentation study found
> coordinate fragmentation is ~1% and does not grow with divergence; a tolerance is a characterized,
> ~30-line fallback if ever needed. Details in
> [`annotate-implementation-notes.md`](./annotate-implementation-notes.md) → "Block-compaction
> open-question resolutions" → issue 1.

**Now:** two per-genome instances cluster only when their whole crossing matches *bit-exactly* —
`(feature_type, [(block_id, cons_start, cons_end, strand_on_consensus), …])` (`ClusterKey`). One indel
near a feature boundary nudges a consensus endpoint by a few bp, which splits one biological gene into
several near-identical clusters, each scored (and possibly dropped) on its own smaller support.

**Why it matters:** this is the single policy choice most likely to make the block view under-report
support on real, divergent genomes. The node-level table still holds everything, so nothing is *lost*
— but the compacted summary can look artificially fragmented.

**Options:**
- **Exact only (status quo).** Predictable and simple; the risk above is real on diverse inputs.
- **Bounded coordinate tolerance.** Treat endpoints within ±δ bp as equal (cluster by rounded/binned
  coordinates, or merge clusters whose termini are within δ). Needs a rule for transitivity/chaining
  and a default δ.

**Recommendation — measure before deciding.** On the klebs / *E. coli* graphs, run a throwaway script
that counts how many biological genes (grouped by GFF `ID`/name) fragment into **more than one**
block-level cluster *purely* because of boundary indels, and how many of those fragments fall below a
0.9 threshold that the whole gene would clear. If that count is small, ship exact-only for v1 and note
the limitation; if it is large, a tolerance is worth the added complexity. This is the concrete
real-data check called for by [`annotate.md`](./annotate.md) §8 ("whether to allow a small coordinate
**tolerance** … is an open refinement").

## 2. Per-segment support pools features that only *look* alike

> **Resolved (2026-07-06): document as intended.** Empirically rare (≤20 placements pool >1 CDS id
> across the study's five datasets) and informational only (gating is per-crossing). Keep the key;
> document that per-segment support is a property of a placement, not a feature.

**Now:** `n_support_segment` is tallied over `segment_support_sets`, keyed on
`(feature_type, block_id, cons_start, cons_end, strand)` with **no feature identity**
(`compact.rs`). Two genuinely different genes of the same type that happen to place a segment at
identical block-consensus coordinates are counted together.

**Why it matters:** it can inflate the reported reproducibility of a segment. Rare in practice (needs
an exact coordinate coincidence between distinct features), but silent when it happens.

**Decision:** either (a) document this as intended — "per-segment support is a property of a
placement, not of a feature" — or (b) add feature identity (e.g. the cluster key, or the base feature
id) to `SegmentKey` so only same-feature segments pool. (a) keeps the "shared body reports the same
support in every crossing" behaviour that motivated the current key; (b) is stricter. Pick one and
state it.

## 3. No provenance / drill-down in the block output

> **Resolved (2026-07-06): defer for v1.** Drill-down means going back to the node-level table; the
> nested-JSON supporters writer lands before the disagreement use case is advertised. P5.3 docs must
> state the node-table drill-down explicitly.

**Now:** a block row reports `M`/`N` and a consensus name/attributes, but not *which* genomes support
or dissent. "Where do genomes disagree about this feature?" is one of the stated use cases, and the
CSV can't answer it.

**Decision:** confirm the deferral is acceptable for v1, and land the nested-JSON writer (already
deferred in [`annotate.md`](./annotate.md) §8) that carries the per-genome supporters behind each
consensus row before the feature is advertised as covering the disagreement use case. Until then, note
in the docs that drill-down means going back to the node-level table.

## 4. Genome identity is a display string, not a `PathId`

> **Resolved (2026-07-06): keyed on `PathId` (implemented).** `LiftedAnnotation` gained a `path_id`
> field and compaction now keys its whole internal pipeline on it (the `genome` string is kept only
> for the node CSV). Removes the collision risk and the per-row `String` clones.

**Now:** the whole compaction keys genomes on a `String` label — `path.name()` falling back to
`path.id().to_string()` — produced identically by the lift and by `genome_label` (so the `M`/`N` join
is correct today).

**Why it matters:** it is a robustness/efficiency point rather than a policy one. Two paths with the
same name (or empty names) would silently merge into one genome, and the label is cloned per row at
genome scale. Keying the internal pipeline on `PathId` (resolving names only at output) removes the
collision risk and the allocation. Low priority; fold in only if the pipeline is touched for another
reason.

## How to close these out

1. Fix the float-rounding off-by-one in `min_count` — **done** (epsilon before `.ceil()`, with a
   regression test).
2. Run the real-data fragmentation measurement (§1) — **done**. The five-species study
   (`tmp/block-annotations/`, note `n00_block_compaction_fragmentation.typ`) decided **exact-only for
   v1**; tolerance stays a ~30-line fallback in one strategy.
3. Make the calls on §2–§4 and record each in
   [`annotate-implementation-notes.md`](./annotate-implementation-notes.md) — **done** (§2 document as
   intended, §3 defer with node-table drill-down, §4 keyed on `PathId`, implemented).

**All four resolved (2026-07-06).** The remaining work is P5.3: write the user docs against this
settled behaviour (exact-only clustering; per-segment support is placement-level; drill-down via the
node-level table).

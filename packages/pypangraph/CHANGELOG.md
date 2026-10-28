## 1.1.1

- Fix duplicated GFA links in `junction_context_gfa`: a link traversed in both orientations within a junction (e.g. accessory content inverted in some isolates) is now written once, and `L` lines are emitted in a deterministic order.

## 1.1.0

Adds **backbone-junction analysis**.

New `pypangraph.junctions` sub-package built around `BackboneJunctions`, which splits each path at core-block boundaries and exposes:

- `stats()`: per-edge summary DataFrame (number of categories, accessory length, occupied isolates, ...).
- `positions()`: genomic coordinates of the flanking core blocks per (edge, isolate).
- `sequences(edge)`: co-oriented `Bio.SeqRecord`s spanning the junction, ready to write to FASTA.

## 1.0.1

Make pypangraph compatible with Python 3.9, by removing `match` syntax.

## 1.0.0

Initial release for pangraph 1.0.0.

For older versions, see https://github.com/mmolari/pypangraph

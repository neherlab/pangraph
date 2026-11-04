# Graph loading and validation

`Pangraph.from_json` turns a pangraph JSON file (optionally gzipped) into an
in-memory graph:

1. **Read.** The file is decompressed when its name ends in `.json.gz` and read
   as bytes.
2. **Parse and validate.** The bytes are validated into the typed model in
   `pypangraph/model.py` with `PangraphData.model_validate_json`. A single call
   parses the JSON and checks structure, field types, the `strand` enum and the
   non-negative integer ranges; a mismatch raises `PangraphLoadError`. There is
   no separate validation step.
3. **Construct.** The typed `PangraphData` is split into the `paths`, `blocks`
   and `nodes` collections that make up a `Pangraph`.

## Why loading parses into typed models

Validation is the dominant cost of loading. On a mid-sized graph (664 blocks,
6817 nodes) a pure-Python JSON-schema pass over the parsed data takes seconds,
because the validator walks every node, edit and position in interpreted Python.

`pydantic` parses and validates JSON in its compiled core, building the typed
models in `model.py` in one call. On the graph above this replaces a multi-second
validated load with about a hundred milliseconds. The models carry the schema's
constraints (`Uint` for non-negative integers, a single-character `alt`, the
`strand` enum), so the accepted and rejected graphs match the schema. A malformed
document (invalid JSON) is reported as a read failure, distinct from a schema
violation.

The models are the internal representation. `Pangraph.from_json` parses into
them; `Pangraph(pan)` also accepts a plain dict, which is validated with
`PangraphData.model_validate`, so building a graph from an in-memory dict enforces
the same schema as loading from a file.

## Reproducing the numbers

`benchmarks/bench_load` measures each phase and every validation engine that is
installed:

```
python3 benchmarks/bench_load                     # tests/data/staph.json.gz
python3 benchmarks/bench_load path/to/graph.json  # a specific graph
```

Install `jsonschema jsonschema-rs msgspec pydantic` to compare all engines on one
machine.

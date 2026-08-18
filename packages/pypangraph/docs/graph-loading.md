# Graph loading and validation

`Pangraph.from_json` turns a pangraph JSON file (optionally gzipped) into an
in-memory graph:

1. **Read.** The file is decompressed when its name ends in `.json.gz` and read
   as bytes.
2. **Decode and validate.** The bytes are decoded straight into the typed model
   in `pypangraph/model.py` with `msgspec`. A single pass checks structure,
   field types, the `strand` enum and the non-negative integer ranges; a mismatch
   raises `PangraphLoadError`. There is no separate validation step.
3. **Construct.** The typed `PangraphData` is split into the `paths`, `blocks`
   and `nodes` collections that make up a `Pangraph`.

## Why loading decodes into typed models

Validation is the dominant cost of loading. On a mid-sized graph (664 blocks,
6817 nodes) a pure-Python JSON-schema pass over the parsed data takes seconds,
because the validator walks every node, edit and position in interpreted Python.

`msgspec` decodes the JSON bytes directly into the compact structs in
`model.py`, so parsing and validation happen together in compiled code and no
intermediate `dict` tree is built. On the graph above this replaces a multi-second
validated load with about twenty milliseconds, faster than parsing into a `dict`
alone. The models carry the schema's constraints (`Uint` for non-negative
integers, a single-character `alt`, the `strand` enum), so the accepted and
rejected graphs match the schema.

The models are the internal representation. `Pangraph.from_json` decodes into
them; `Pangraph(pan)` also accepts a plain dict, which is validated and converted
with `msgspec.convert`, so building a graph from an in-memory dict enforces the
same schema as loading from a file.

## Reproducing the numbers

`benchmarks/bench_load` measures each phase and every validation engine that is
installed:

```
python3 benchmarks/bench_load                     # tests/data/staph.json.gz
python3 benchmarks/bench_load path/to/graph.json  # a specific graph
```

Install `jsonschema jsonschema-rs msgspec pydantic` to compare all engines on one
machine.

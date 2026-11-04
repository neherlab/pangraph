# Graph loading and validation

`Pangraph.from_json` turns a pangraph JSON file (optionally gzipped) into an
in-memory graph in three steps:

1. **Read and parse.** The file is decompressed when its name ends in `.json.gz`
   and parsed with the standard library `json` module.
2. **Validate.** The parsed object is checked against the JSON schema in
   `pypangraph/pangraph_schema.py`, which is generated from the Rust types that
   produce the file. Validation failures raise `PangraphLoadError` with the
   offending constraint.
3. **Construct.** The validated object is split into the `paths`, `blocks` and
   `nodes` collections that make up a `Pangraph`.

## Why validation uses jsonschema-rs

Validation is the dominant cost of loading. On a mid-sized graph (664 blocks,
6817 nodes) parsing takes tens of milliseconds while a pure-Python JSON-schema
pass over the same data takes seconds: the validator walks every node, edit and
position in interpreted Python.

The loader validates with [`jsonschema-rs`](https://pypi.org/project/jsonschema-rs/),
a Rust-backed validator for the same schema. It compiles the schema once at import
time (`_VALIDATOR` in `pypangraph/class_graph.py`) and reuses the compiled
validator for every load. On the graph above this reduces the validation step from
seconds to about ten milliseconds, so the load becomes bounded by decompression and
parsing rather than validation. The accepted and rejected graphs are unchanged:
`jsonschema-rs` enforces the same schema, including required fields, value types,
the `strand` enum, and non-negative integers.

## Reproducing the numbers

`benchmarks/bench_load` measures each phase and every validation engine that is
installed:

```
python3 benchmarks/bench_load                     # tests/data/staph.json.gz
python3 benchmarks/bench_load path/to/graph.json  # a specific graph
```

Install `jsonschema jsonschema-rs msgspec pydantic` to compare all engines on one
machine.

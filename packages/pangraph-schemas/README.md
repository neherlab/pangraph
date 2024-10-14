# Pangraph schemas

This directory contains [JSON Schema](https://json-schema.org/) definitions for Pangraph types in the current branch of the project. They are automatically generated from Rust types using [`/packages/pangraph-cli/src/build.rs`](../pangraph-cli/src/build.rs) build script (See: [Cargo build script docs](https://doc.rust-lang.org/cargo/reference/build-scripts.html)), which runs every time Pangraph CLI is built.

Note that these definitions correspond to the current commit in the Pangraph source code repo. As project changes, the type definitions can also change and not necessarily reflect the format used in the particular version of Pangraph you are using.

In order to generate definitions for the version of pangraph you are currently using, you can run

```bash
pangraph schema -o Pangraph.schema.json
```

or

```bash
pangraph schema -o Pangraph.schema.yml
```

You can specify either JSON or YAML for the output file format.

From there, you can use the generated JSON schema as you'd like. For example, you can use `datamodel-code-generator` to generate corresponding Python classes and use them in your programs:

```bash
cd pangraph-schemas/

pip3 install dacite datamodel-code-generator pydantic

# Generate Python classes
datamodel-codegen --input-file-type "jsonschema" --input "Pangraph.schema.json" --output-model-type "dataclasses.dataclass" --enum-field-as-literal=all --output "python/Pangraph.py"

# Run the example using the generated Python classes
python3 python/example.py path/to/your/pangraph/output.json

```

In this example the generated `Pangraph.py` file will contain the Python dataclasses. The example reads the Pangraph CLI file and casts it to the generated `Pangraph` dataclass types using `dacite` library. You can then access to the data in a type-safe manner, and most text editors should also provide code completions.

JSON Schema is a popular type definition format, so most programming languages have tools and libraries to work with it. Feel free explore some of the libraries in your favourite language! There are also multi-language tools and even online generators like [Qucktype](https://quicktype.io/). Quicktype also has a Node.js CLI:

```bash
npx -y quicktype --src-lang schema --src "Pangraph.schema.json" --lang python --python-version 3.7 --just-types --top-level "_PangraphSchemaRoot" --out "python/Pangraph.py"
```

> ⚠️ Note the due to a number of existing programming languages, large variety of tools, as well as differences in the code they may generate, the dev team cannot guarantee correctness of the generated code, and neither cannot provide technical support for any of the third-party tools.
> 
> Think of the provided JSON Schema definitions as a helper utility and a machine-readable documentation for the file formats. It is a starting point rather than a complete solution. Depending on the language, quality of your code generator and your goals, you might need to experiment a little to make things work.
> 
> That being said, if you want to suggest an improvement to how the JSON Schemas themselves are generated, feel free to open a GitHub issue.

SHELLFLAGS := -euo pipefail
SHELL := bash
.ONESHELL:

SCHEMA_FILE := packages/pangraph-schemas/Pangraph.schema.json
SCHEMA_FILE_PYPANGRAPH := packages/pypangraph/pypangraph/Pangraph.schema.json

OUTPUT_FILE_PYPANGRAPH := packages/pypangraph/pypangraph/Pangraph_model.py
OUTPUT_FILE_EXAMPLE := packages/pangraph-schemas/python/Pangraph_model.py

.PHONY: requirements clean schema codegen

all: clean schema codegen

# Install required dependencies
requirements:
	pip3 install -qr requirements.txt

# Remove generated files
clean:
	rm -rf $(SCHEMA_FILE) $(SCHEMA_FILE_PYPANGRAPH) $(OUTPUT_FILE_PYPANGRAPH) $(OUTPUT_FILE_EXAMPLE)

# Generate JSON schema using Pangraph CLI
schema: $(SCHEMA_FILE) $(SCHEMA_FILE_PYPANGRAPH)

$(SCHEMA_FILE):
	cargo run -q --bin=pangraph -- schema -o "$(SCHEMA_FILE)"

# Generate Python code from JSON schema
codegen: $(OUTPUT_FILE_PYPANGRAPH) $(OUTPUT_FILE_EXAMPLE)

$(OUTPUT_FILE_EXAMPLE): $(SCHEMA_FILE)
	datamodel-codegen \
		--input-file-type "jsonschema" \
		--input "$(SCHEMA_FILE)" \
		--output-model-type "dataclasses.dataclass" \
		--enum-field-as-literal=all \
		--output "$(OUTPUT_FILE_EXAMPLE)"

$(SCHEMA_FILE_PYPANGRAPH): $(SCHEMA_FILE)
	cp $(SCHEMA_FILE) $(SCHEMA_FILE_PYPANGRAPH)

$(OUTPUT_FILE_PYPANGRAPH): $(OUTPUT_FILE_EXAMPLE)
	cp $(OUTPUT_FILE_EXAMPLE) $(OUTPUT_FILE_PYPANGRAPH)

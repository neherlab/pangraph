SHELLFLAGS := -euxo pipefail
SHELL := bash
.ONESHELL:

SCHEMA_FILE="packages/pangraph-schemas/Pangraph.schema.json"
OUTPUT_FILE_PYPANGRAPH="packages/pypangraph/pypangraph/Pangraph_model.py"
OUTPUT_FILE_EXAMPLE="packages/pangraph-schemas/python/Pangraph.py"

.PHONY: requirements clean schema codegen

all: clean schema codegen

# Install required dependencies
requirements:
	@set -euxo pipefail
	@pip3 install -qr requirements.txt

# Remove generated files
clean:
	@set -euxo pipefail
	@rm -rf $(OUTPUT_FILE_PYPANGRAPH) $(OUTPUT_FILE_EXAMPLE)

# Generate JSON schema using Pangraph CLI
schema: $(SCHEMA_FILE)

$(SCHEMA_FILE):
	@set -euxo pipefail
	@cargo run -q --bin=pangraph -- schema -o "$(SCHEMA_FILE)"

# Generate Python code from JSON schema
codegen: $(OUTPUT_FILE_PYPANGRAPH) $(OUTPUT_FILE_EXAMPLE)

$(OUTPUT_FILE_PYPANGRAPH): $(SCHEMA_FILE)
	@set -euxo pipefail
	@datamodel-codegen \
		--input "$(SCHEMA_FILE)" \
		--output "$(OUTPUT_FILE_PYPANGRAPH)" \
		--input-file-type "jsonschema" \
		--output-model-type "dataclasses.dataclass" \
		--enum-field-as-literal=all \

$(OUTPUT_FILE_EXAMPLE): $(OUTPUT_FILE_PYPANGRAPH)
	@cp $(OUTPUT_FILE_PYPANGRAPH) $(OUTPUT_FILE_EXAMPLE)

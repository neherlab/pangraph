#!/bin/bash

set -euxo pipefail

python3 tools/release_checks.py --release $1 --changelog "CHANGELOG.md" --project "Project.toml"
python3 tools/extract_release_notes.py "CHANGELOG.md"
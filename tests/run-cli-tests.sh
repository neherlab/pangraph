#!/usr/bin/env bash

set -euxo pipefail

: "${1:? ${0}: Path to pangraph executable is required}"

function abspath() {
  readlink -m "$1"
}

# shellcheck disable=SC2155
export THIS_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

# shellcheck disable=SC2155
export PROJECT_ROOT_DIR="$(abspath "${THIS_DIR}/..")"

PANGRAPH="${1}"
IN_DIR="${PROJECT_ROOT_DIR}/tests/in"
OUT_DIR="${PROJECT_ROOT_DIR}/tests/out"

mkdir -p "$OUT_DIR"

# Test version
${PANGRAPH} --version

# Test long help
${PANGRAPH} --help
${PANGRAPH} build --help
${PANGRAPH} export --help
${PANGRAPH} generate --help
${PANGRAPH} marginalize --help
${PANGRAPH} polish --help

# Test short help
${PANGRAPH} -h
${PANGRAPH} build -h
${PANGRAPH} export -h
${PANGRAPH} generate -h
${PANGRAPH} marginalize -h
${PANGRAPH} polish -h

# Test pangraph generate
${PANGRAPH} generate -v -s 100 -r 1e-1 -m 1e-3 -d 5e-2 -i 1e-2 -t 5 "$IN_DIR/randseqs.fa" >"$OUT_DIR/input.fa"

# Test pangraph build - minimap asm20 no energy
${PANGRAPH} build -v -c -k minimap2 -s 20 -a 0 -b 0 "$OUT_DIR/input.fa" >"$OUT_DIR/test1.json"

# Test pangraph build - mash
${PANGRAPH} build -v -c -d mash "$OUT_DIR/input.fa" >"$OUT_DIR/test2.json"

# Test pangraph build - mmseqs
${PANGRAPH} build -v -c -k mmseqs -K 8 "$OUT_DIR/input.fa" >"$OUT_DIR/test3.json"

# "Test pangraph polish
${PANGRAPH} polish -v -c -l 10000 "$OUT_DIR/test1.json" >"$OUT_DIR/polished.json"

# Test pangraph GFA export
${PANGRAPH} export -v -o "$OUT_DIR/export" "$OUT_DIR/test1.json"

# Test pangraph PanX export
${PANGRAPH} export -v --ng --px -o "$OUT_DIR/export" "$OUT_DIR/test1.json"

# Test pangraph marginalize
${PANGRAPH} marginalize -v -o "$OUT_DIR/marginalize" "$OUT_DIR/test1.json"

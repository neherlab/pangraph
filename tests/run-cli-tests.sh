#!/bin/bash
set -euxo pipefail

# test that mash, mafft, mmseqs and fasttree are available in path
echo "mash version:"
mash --version
echo "mafft version:"
mafft --version
echo "mmseqs version:"
mmseqs --help | grep "Version"
echo "fasttree help:"
fasttree -help

# test pangraph commands help
pangraph help build
pangraph help polish
pangraph help export
pangraph help marginalize

# create input data
TESTDIR="tests/data"
mkdir -p "$TESTDIR"

julia --project=. -e "using PanGraph; \
    sequences=[PanGraph.Simulation.randseq(30000) for _ in 1:6]; \
    open(\"$TESTDIR/randseqs.fa\", \"w\") do io ;\
        for (i,sequence) in enumerate(sequences) ; \
            PanGraph.Simulation.write_fasta(io, \"isolate_$(i)\", sequence) ; \
        end ; \
    end"

# test pangraph commands

echo "Test pangraph help"
pangraph help
pangraph help help
pangraph help build
pangraph help generate
pangraph help version
pangraph help export
pangraph help marginalize
pangraph help polish

echo "Test pangraph version"
pangraph version

echo "Test pangraph generate"
pangraph generate -s 100 -r 1e-1 -m 1e-3 -d 5e-2 -i 1e-2 -t 5 "$TESTDIR/randseqs.fa" > "$TESTDIR/input.fa"

echo "Test pangraph build - minimap asm20 no energy"
export JULIA_NUM_THREADS=4
pangraph build -v -c -k minimap2 -s 20 -a 0 -b 0 "$TESTDIR/input.fa" > "$TESTDIR/test1.json"

echo "Test pangraph build - mash"
pangraph build -v -c -d mash "$TESTDIR/input.fa" > "$TESTDIR/test2.json"

echo "Test pangraph build - mmseqs"
export JULIA_NUM_THREADS=1
pangraph build -v -c -k mmseqs -K 8 "$TESTDIR/input.fa" > "$TESTDIR/test3.json"

echo "Test pangraph polish"
pangraph polish -c -l 10000 "$TESTDIR/test1.json" > "$TESTDIR/polished.json"

echo "Test pangraph GFA export"
pangraph export -o "$TESTDIR/export" "$TESTDIR/test1.json"

echo "Test pangraph PanX export"
pangraph export -ng -pX -o "$TESTDIR/export" "$TESTDIR/test1.json"

echo "Test pangraph marginalize"
pangraph marginalize -v -o "$TESTDIR/marginalize" "$TESTDIR/test1.json"

# remove generated files
rm -r $TESTDIR

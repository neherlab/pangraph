# This script convers a genbank file to a fasta file. Only the record with the longes
# sequence is saved. The fasta record id is the same as the genbank filename.
# The script checks that the original sequence id is compatible with this name.

from Bio import SeqIO
import argparse
import re
import numpy as np

def accnum(fname):
    """Extract the .gbk file name without extension"""
    p = re.compile("/?([^/]+)\.gbk$")
    return re.search(p, fname).groups()[0]

if __name__ == "__main__":

    # parse arguments
    parser = argparse.ArgumentParser(description="convert a genbank file to a fasta file. Only the record with the longes sequence is saved.")
    parser.add_argument("--gbk", help="input genbank file.", type=str)
    parser.add_argument("--fa", help="output fasta file.", type=str)
    args = parser.parse_args()

    # read reacords in the input genbank file
    records = list(SeqIO.parse(args.gbk, "genbank"))

    # pick the longest record
    idx = 0
    NR = len(records)
    if NR > 1:
        print(f"{NR} records found. Picking the longest.")
        Ls = [len(r.seq) for r in records]
        idx = np.argmax(Ls)
        print(f"Picked record {records[idx].id} with length {Ls[idx]//1000} kbp")
    r = records[idx]

    # check that id of the read is compatible with filename
    acc = accnum(args.gbk)
    assert r.id.startswith(acc), f"read id {r.id} and filename {acc} are incompatible."

    # use only name of the file as id of the read
    r.id = acc
    r.name = ""
    r.description = ""
    
    # write the fasta file
    SeqIO.write([r], args.fa, "fasta")
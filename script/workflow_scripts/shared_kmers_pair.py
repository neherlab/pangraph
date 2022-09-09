import argparse
import json
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Evaluates the fraction of kmers shared by two sequences
        and returns it in a """
    )
    parser.add_argument("--s1", help="input fasta file for first isolate", type=str)
    parser.add_argument("--s2", help="input fasta file for second isolate", type=str)
    parser.add_argument("--k", help="length of kmers", type=int)
    parser.add_argument("--json", help="output json file", type=str)
    return parser.parse_args()


def kmers_set(biopython_seq, kmer_size):
    """Given a biopython seq object and k-mer size returns the set of kmers"""
    seq = str(biopython_seq.seq)
    L = len(seq)
    kmers = [seq[i : i + kmer_size] for i in range(L - kmer_size)]
    return set(kmers)


def pack_results(s, pair, K, KS):
    return {
        "id": s.id,
        "pair": pair,
        "ktot": K,
        "shared": KS,
        "f": KS / K,
    }


if __name__ == "__main__":

    args = parse_args()

    # load sequences
    s1 = SeqIO.read(args.s1, format="fasta")
    s2 = SeqIO.read(args.s2, format="fasta")

    # evaluate kmer sets
    k1 = kmers_set(s1, kmer_size=args.k)
    k2 = kmers_set(s2, kmer_size=args.k)

    # evaluate shared kmers
    shared_k = k1.intersection(k2)

    # evaluate kmer numbers
    K1, K2, KS = len(k1), len(k2), len(shared_k)

    # pack and save results
    pair = f"{s1.id}-{s2.id}"
    r1 = pack_results(s1, pair, K1, KS)
    r2 = pack_results(s2, pair, K2, KS)
    with open(args.json, "w") as f:
        json.dump([r1, r2], f)

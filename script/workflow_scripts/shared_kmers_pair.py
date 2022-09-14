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


def pack_results(s, pair, K, KS, KT):
    return {
        "id": s.id, # strain id
        "pair": pair, # pair id ("s1-s2")
        "n_strain": K, # number of strain k-mers
        "n_union": KT, # number of kmer union
        "n_intersection": KS, # number of k-mer intersection
        "f": KS / K, # fraction of shared k-mers on the single strain
        "jaccardi": KS / KT, # jaccardi index
    }


if __name__ == "__main__":

    args = parse_args()

    # load sequences
    s1 = SeqIO.read(args.s1, format="fasta")
    s2 = SeqIO.read(args.s2, format="fasta")

    # evaluate kmer sets
    k1 = kmers_set(s1, kmer_size=args.k)
    k2 = kmers_set(s2, kmer_size=args.k)

    # evaluate number of shared and total kmers
    K1, K2 = len(k1), len(k2)
    KS = len(k1.intersection(k2))
    KT = len(k1.union(k2))

    # pack and save results
    pair = f"{s1.id}-{s2.id}"
    r1 = pack_results(s1, pair, K1, KS, KT)
    r2 = pack_results(s2, pair, K2, KS, KT)
    with open(args.json, "w") as f:
        json.dump([r1, r2], f)

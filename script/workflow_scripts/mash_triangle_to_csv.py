import argparse
import re
import numpy as np
import pandas as pd


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="""Turns the results of mash-triangle into a csv file with 
        pairs of strains as index and mash distances as values."""
    )
    parser.add_argument("--mash_tri", help="output of mash triangle", type=str)
    parser.add_argument("--csv", help="output csv file", type=str)
    return parser.parse_args()


def parse_mash_triangle(fname: str):

    # parse file
    with open(fname, "r") as f:
        lines = f.readlines()

    # n. isolates
    N = int(lines[0].strip())

    # distance matrix and strain names
    M = np.zeros((N, N))
    strains = []

    # capture strain names and distances
    for i, l in enumerate(lines[1:]):
        vals = l.split()
        name = re.search(r"/([^/]+)\.fa$", vals[0]).group(1)
        strains.append(name)
        for j, v in enumerate(vals[1:]):
            M[i, j] = float(v)
            M[j, i] = float(v)

    return strains, M


def distance_mat_to_df(S, M):

    # capture pairs of strains and their distance
    strains = []
    distances = []
    for i, s1 in enumerate(S):
        for j, s2 in enumerate(S):
            strains.append((s1, s2))
            distances.append(M[i, j])
    
    # save in a dataframe
    idx = pd.MultiIndex.from_tuples(strains, names=["strain_1", "strain_2"])
    df = pd.DataFrame({"mash_dist": distances}, index=idx)

    return df


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # parse mash distance file
    strains, M = parse_mash_triangle(args.mash_tri)

    # turn into multi-index dataframe
    df = distance_mat_to_df(strains, M)

    # save distance dataframe
    df.to_csv(args.csv)

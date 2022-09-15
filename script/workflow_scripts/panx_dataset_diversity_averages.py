import argparse
import re
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Parses all species diversity dataframes, evaluates average
        statistics and returns them in a csv dataframe."""
    )
    parser.add_argument(
        "--dfs", help="input species dataframes (csv)", nargs="+", type=str
    )
    parser.add_argument("--out", help="output csv file", type=str)
    return parser.parse_args()


def species_statistics(df):
    """Evaluate summary statistics for a species dataframe"""

    # core genes dataframe
    mask = df["core"] & df["non_cons_fraction"].notna()
    cdf = df[mask]

    # function to evaluate weighed averages
    def w_avg(dflab, wlab):
        return np.sum(cdf[dflab] * cdf[wlab]) / np.sum(cdf[wlab])

    return {
        "avg_core_noncons_fraction": w_avg("non_cons_fraction", "len_nogap"),
        "avg_core_pairwise_divergence": w_avg("pairwise_divergence", "len_nogap"),
        "avg_core_shared_kmer_fraction": w_avg("avg_kmer_shared_fract", "len_avg"),
        "avg_core_jaccardi_index": w_avg("avg_kmer_jaccardi", "len_avg"),
        "pangenome_len": df["len_avg"].sum(),
        "core_genome_len": cdf["len_avg"].sum(),
        "n. genes": len(df),
        "n. core genes": len(cdf),
    }


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # evaluate averages and create summary dataframe
    avg_df = []
    for df_file in args.dfs:

        # load dataframe
        df = pd.read_csv(df_file)

        # evaluate summary statistics
        info = species_statistics(df)

        # parse species name from filename
        info["species"] = re.search(r"\/([^\/]+)\.csv$", df_file).group(1)

        avg_df.append(info)

    # convert to dataframe
    avg_df = pd.DataFrame(avg_df).set_index("species")

    # save dataframe
    avg_df.to_csv(args.out)

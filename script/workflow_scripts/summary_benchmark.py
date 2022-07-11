import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# patterns to capture from the txt file
patterns = {
    "species" : "species = ([\S]+)", # string
    "kind" : "kind = ([\S]+)", # string
    "n-isolates" : 'Command being timed: "(.*)"',
    "mem" : "Maximum resident set size \(kbytes\): ([\d]+)", # int
    "wall-time" : "Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\S]+)", # time
    "usr-time" : "User time \(seconds\): ([\S]+)", # float
    "sys-time" : "System time \(seconds\): ([\S]+)", # float
    "cpu-percent" : "Percent of CPU this job got: ([\d]+)%", # int
    "command" : 'Command being timed: "(.*)"', #str
}

# time in format hh:mm:ss.ss into seconds
def parsetime(x):
    s = [float(k) for k in x.split(":")]
    t = s[-1] + 60 * s[-2]
    if len(s) > 2:
        t += s[-3] * 3600
    return t


# conversion functions to transform captured patterns into values
convert = {
    "species" : str,
    "kind" : str,
    "n-isolates" : lambda x : len([y for y in x.split(" ") if y.startswith("panx_data/")]),
    "mem" : lambda x : int(x) / 1e6,
    "wall-time" : parsetime,
    "usr-time" : lambda x: float(x) / (3600),
    "sys-time" : lambda x: float(x) / (3600),
    "cpu-percent" : int,
    "command" : lambda x : " ".join([y for y in x.split(" ") if not y.startswith("panx_data/")]),
}

# given a pattern extract the first occurrence from the file
def extract_pattern(m, lines):
    x = re.findall(re.compile(m), lines)
    if len(x) == 0:
        return ""
    else:
        return x[0]

# extract all patterns from a list of files
def extract_info(f_name):
    # read file in a single string
    with open(f_name, "r") as f:
        lines = f.read()
    # capture all patterns
    captures = {}
    for k, m in patterns.items():
        val = extract_pattern(m, lines)
        captures[k] = convert[k](val)
    return captures

# plot with benchmark statistics
def plot_benchmark(df):

    # modify a copy of the dataframe
    df = df.copy()

    # from full species name to abbreviation 
    def short_name(x):
        y = x.split("_")
        y[0] = y[0][0]
        y = [str.capitalize(s) for s in y]
        return ". ".join(y)
    df["species"] = df["species"].apply(short_name)

    # order the species
    order = df.groupby("species")["n-isolates"].mean().sort_values().index.to_numpy()
    df["species"] = pd.Categorical(df["species"], categories=order, ordered=True)

    # create figure
    N = 4
    fig, axs = plt.subplots(N,1, figsize=(8,N*4))

    # plot n. isolates
    ax = axs[0]
    sns.barplot(x="species", y="n-isolates", data=df, ax=ax)

    # plot wall-time
    ax=axs[1]
    df["wall-time (h)"] = df["wall-time"] / 3600
    sns.barplot(x="kind", hue="species", y="wall-time (h)", data=df, ax=ax)
    ax.set_yscale("log")

    # plot max resident size
    ax = axs[2]
    sns.barplot(x="kind", hue="species", y="mem", data=df, ax=ax)
    ax.set_ylabel("max resident size (Gb)")

    # plot cpu percent
    ax = axs[3]
    sns.barplot(x="kind", hue="species", y="cpu-percent", data=df, ax=ax)


    return fig, axs



if __name__ == "__main__":

    # parse arguments
    parser = argparse.ArgumentParser(description="creates a table with summary statistics from the benchmark data")
    parser.add_argument("out_csv", help="name of the output csv table", type=str)
    parser.add_argument("out_pdf", help="name of the output plots", type=str)
    parser.add_argument("benchmarks", nargs="*", help="list of input summary txt files", type=str)
    args = parser.parse_args()

    # extract info from input files
    in_files = list(args.benchmarks)
    infos = [extract_info(fn) for fn in in_files]

    # create and export dataframe
    df = pd.DataFrame(infos)
    df.to_csv(args.out_csv, index=False)

    # plot results
    fig, axs = plot_benchmark(df)
    plt.tight_layout()
    plt.savefig(args.out_pdf)
    plt.close(fig)

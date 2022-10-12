# Script to parse the output size benchmark txt files and collect
# the results in a csv dataframe.

import argparse
import re
import pandas as pd


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="""parses the output size benchmark txt files and collect
        the results in a csv dataframe."""
    )
    parser.add_argument("--txt", nargs="+", help="input txt benchmark files", type=str)
    parser.add_argument("--csv", help="output csv dataframe", type=str)
    return parser.parse_args()


# patterns to capture from the txt file
patterns = {
    "n-isolates": "N = ([\d]+)",  # number
    "length": "L = ([\d]+)",  # number
    "trial": "trial = ([\d]+)",  # number
    "mem": "Maximum resident set size \(kbytes\): ([\d]+)",  # int
    "wall-time": "Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\S]+)",  # time
    "cpu-percent": "Percent of CPU this job got: ([\d]+)%",  # int
    "command": 'Command being timed: "(.*)"',  # str
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
    "n-isolates": int,
    "length": int,
    "trial": int,
    "mem": lambda x: int(x) / 1e6,
    "wall-time": parsetime,
    "cpu-percent": int,
    "command": str,
}

# given a pattern extract the first occurrence from a string
def extract_pattern(m, content):
    x = re.findall(re.compile(m), content)
    if len(x) == 0:
        return ""
    else:
        return x[0]


# given a text file it reads it and captures all of the
def parse_benchmark(fname):

    # read all the file in a single string
    with open(fname, "r") as f:
        content = f.read()

    # capture all patterns in a dictionary
    captures = {}
    for k, m in patterns.items():
        val = extract_pattern(m, content)
        captures[k] = convert[k](val)

    return captures


if __name__ == "__main__":

    # parse args
    args = parse_args()

    # parse text files
    df = []
    for fname in args.txt:
        df.append(parse_benchmark(fname))

    # build and save dataframe
    df = pd.DataFrame(df)
    df.to_csv(args.csv, index=False)

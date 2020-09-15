#!/usr/bin/env python3
"""
script to process our end repair log files for plotting
"""
import ast
import argparse

from io          import StringIO
from Bio         import AlignIO
from collections import defaultdict

import numpy as np
import matplotlib.pylab as plt

level_preset = "+++LEVEL="
level_offset = len(level_preset)

score_preset = "SCORE="
score_offset = len(score_preset)

msaln_preset = "ALIGNMENT="
msaln_offset = len(msaln_preset)

def unpack(line):
    offset = [line.find(";")]
    offset.append(line.find(";", offset[0]+1))
    offset.append(line.find(";", offset[1]+1))

    align = StringIO(ast.literal_eval(line[msaln_offset+1+offset[0]:offset[1]]).decode('utf-8'))
    align = next(AlignIO.parse(align, 'fasta'))
    score = float(line[offset[1]+1+score_offset:offset[2]])

    return align, score

# TODO: remove hardcoded numbers
def save_aln_examples(results):
    stats = results[(500,1000)]
    nums  = [0, 0, 0]
    for score, aln in stats[1]['hits']:
        if len(aln) < 10:
            continue

        if score < 1e-4:
            with open(f"scratch/1/eg_{nums[0]}.fna", 'w+') as fd:
                AlignIO.write(aln, fd, "fasta")
                nums[0] += 1
        elif score < 1e-2:
            with open(f"scratch/2/eg_{nums[1]}.fna", 'w+') as fd:
                AlignIO.write(aln, fd, "fasta")
                nums[1] += 1
        else:
            with open(f"scratch/3/eg_{nums[2]}.fna", 'w+') as fd:
                AlignIO.write(aln, fd, "fasta")
                nums[2] += 1

def main(args):
    results = {}
    for log_path in args:
        stats = defaultdict(lambda: {'hits':[], 'miss': 0})
        level = -1
        with open(log_path) as log:
            for line in log:
                line.rstrip('\n')
                if line[0] == "+":
                    assert line.startswith(level_preset), "check syntax in log file"
                    level = int(line[level_offset:line.find("+++", level_offset)])
                    continue
                if line[0] == ">":
                    if line[1:].startswith("NO MATCH"):
                        stats[level]['miss'] += 1
                        continue
                    if line[1:].startswith("LEN="):
                        aln, score = unpack(line)
                        stats[level]['hits'].append((score, aln))
                        continue
                raise ValueError(f"invalid syntax: {line[1:]}")
        if len(stats) > 0:
            path = log_path.replace(".log", "").split("-")
            e, w = int(path[1][1:]), int(path[2][1:])
            results[(e,w)] = dict(stats)
    return results

def plot(results):
    data   = results[(500,1000)]
    scores = [[] for _ in data.keys()]
    fracs  = np.zeros(len(data.keys()))
    for i, lvl in enumerate(sorted(data.keys())):
        scores[i] = [elt[0] for elt in data[lvl]['hits']]
        fracs[i]  = len(scores[i])/(len(scores[i]) + data[lvl]['miss'])

    cmap   = plt.cm.get_cmap('plasma')
    colors = [cmap(x) for x in np.linspace(0,1,len(scores[:-5]))]
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(fracs[:-5])
    ax1.set_xlabel("tree level")
    ax1.set_ylabel("fraction of good edges")
    for i, score in enumerate(scores[:-5]):
        ax2.plot(sorted(np.exp(score)), np.linspace(0,1,len(score)),color=colors[i],label=f"level={i+1}")
    ax2.set_xlabel("average column entropy")
    ax2.set_ylabel("CDF")
    ax2.legend()

parser = argparse.ArgumentParser(description='process our data log files on end repair')
parser.add_argument('files', type=str, nargs='+')

# -----------------------------------------------------------------------------
# main point of entry

if __name__ == "__main__":
    args = parser.parse_args()
    results = main(args.files)
    save_aln_examples(results)
    plot(results)

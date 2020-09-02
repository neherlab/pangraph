#!/usr/bin/env python3
"""
script to process our end repair log files for plotting
"""
import argparse
from collections import defaultdict

level_preset = "+++LEVEL="
level_offset = len(level_preset)

score_preset = "SCORE="
score_offset = len(score_preset)

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
                        offset = [line.find(";")]
                        offset.append(line.find(";", offset[0]+1))
                        offset.append(line.find(";", offset[1]+1))

                        score    = [None, None]
                        score[0] = float(line[offset[0]+1+score_offset:offset[1]])
                        score[1] = float(line[offset[1]+1+score_offset:offset[2]])
                        stats[level]['hits'].extend(score)
                        continue
                raise ValueError(f"invalid syntax: {line[1:]}")
        if len(stats) > 0:
            path = log_path.replace(".log", "").split("-")
            e, w = int(path[1][1:]), int(path[2][1:])
            results[(e,w)] = dict(stats)
    return results

parser = argparse.ArgumentParser(description='process our data log files on end repair')
parser.add_argument('files', type=str, nargs='+')

if __name__ == "__main__":
    args = parser.parse_args()
    results = main(args.files)

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
    for log_path in args:
        with open(log_path) as log:
            level = -1
            stats = defaultdict(lambda: {'hits':[], 'miss': 0})
            for line in log:
                line.rstrip('\n')
                if line[0] == "+":
                    assert line.startwith(level_preset), "check syntax in log file"
                    level = int(line[level_offset:line.find("+++", level_offset)])
                    continue
                if line[0] == ">":
                    if line[1:] == "NO MATCH":
                        stats[level]['miss'] += 1
                        continue
                    if line[1:].startswith("LEN="):
                        offset = [line.find(";")]
                        offset.append(line.find(";", offset[0]))
                        score[0] = int(line[line.find(score_preset, offset[0])+1+score_offset:offset[1])
                        score[1] = int(line[line.find(score_preset, offset[1])+1+score_offset:)
                        stats[level]['hits'].extend(score)

                raise ValueError(f"invalid syntax: {line}")
        print(stats)

parser = argparse.ArgumentParser(description='process our data log files on end repair')
parser.add_argument('files', type=str, nargs='+')

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.files)

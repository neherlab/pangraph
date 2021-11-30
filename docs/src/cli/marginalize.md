# Marginalize

## Description
Compute all pairwise marginalizations of a multiple sequence alignment pangraph.

## Options
Name | Type | Short Flag | Long Flag | Description
:-------------- | :------- | :------ | :------- | :-------------------------
Output path | String | o | output-path | Path to direcotry where output files will be stored
Reduce paralogs | Boolean | r | reduce-paralog | Collapses coparallel paths through duplicated blocks.

## Arguments
Zero or one pangraph file which must be formatted as a JSON.
If no file path is given, reads from _stdin_.
In either case, the stream can be optionally gzipped.

## Output
Outputs all pairwise graphs to the directory at the user-supplied path.

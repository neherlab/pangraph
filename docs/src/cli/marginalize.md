# Marginalize

## Description
Compute all pairwise marginalizations of a multiple sequence alignment pangraph.

## Options
Name | Type | Short Flag | Long Flag | Description
:-------------- | :------- | :------ | :------- | :-------------------------
Output path | String | o | output-path | Path to direcotry where the output of all pairwise mariginalizations will be stored if supplied
Reduce paralogs | Boolean | r | reduce-paralog | Collapses coparallel paths through duplicated blocks.
Projection strains | String | s | Strains | Collapses the graph structure to only blocks and edges contained by the paths of the supplied strain names. comma seperated, no spaces

## Arguments
Zero or one pangraph file which must be formatted as a JSON.
If no file path is given, reads from _stdin_.
In either case, the stream can be optionally gzipped.

## Output
Outputs all pairwise graphs to the directory at the user-supplied path.

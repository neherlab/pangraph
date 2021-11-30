# Generate

## Description
Generate a simulated multiple sequence alignment pangraph.

## Options
Name | Type | Short Flag | Long Flag | Description
:-------------- | :------- | :------ | :------- | :-------------------------
Mutation rate | Float | m | snp-rate | Rate of mutations per site per genome per generation
HGT rate | Float | r | hgt-rate | Rate of horizontal transfer events per genome per generation
Deletion rate | Float | d | delete-rate | Rate of deletion events per genome per generation
Inversion rate | Float | i | invert-rate | Rate of inversion events per genome per generation
Graph output | String | o | output-path | Path to location to store simulated pangraph
Time | Integer | t | time | Number of generations to simulate before computing sequences and graph

## Arguments
Zero or one fasta file to treat as ancestral sequences.
If no file path is given, reads from _stdin_.
In either case, the stream can be optionally gzipped.
The number and length of sequences determine the population size.

## Output
Outputs all resultant sequences to standard out.
Optionally output the resultant pangraph if a path is supplied by the user.

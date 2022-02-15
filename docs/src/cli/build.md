# Build

## Description
Build a multiple sequence alignment pangraph.

## Options
Name | Type | Short Flag | Long Flag | Description
:-------------- | :------- | :------ | :------- | :-------------------------
minimum length | Integer | l | len | minimum block size for alignment graph (in nucleotides)
block junction cost | Float | a | alpha| energy cost for interblock diversity due to alignment merger,
circular genomes | Boolean | c | circular | toggle if input genomes are circular
pairwise sensitivity | String | s | sensitivity | controls the pairwise genome alignment sensitivity. currently only accepts "5", "10" or "20"
maximum self-maps | Integer | x | max-self-map | maximum number of iterations to perform block self maps per pairwise graph merger
enforce uppercase | Boolean | u | upper-case | toggle if input genomes are set to uppercase characters
distance calculator | String | d | distance-backend | only accepts "native" or "mash"

## Arguments
Expects one or more fasta files.
Multiple records within one file are treated as separate genomes
Fasta files can be optionally gzipped.

## Output
Prints the constructed pangraph as a JSON to _stdout_.

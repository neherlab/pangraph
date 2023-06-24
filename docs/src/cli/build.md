# Build

## Description
Build a multiple sequence alignment pangraph.

## Options
| Name                 | Type    | Short Flag | Long Flag        | Description                                                                                                   |
| :------------------- | :------ | :--------- | :--------------- | :------------------------------------------------------------------------------------------------------------ |
| minimum length       | Integer | l          | len              | minimum block size for alignment graph (in nucleotides)                                                       |
| block junction cost  | Float   | a          | alpha            | energy cost for introducing block partitions due to alignment merger                                          |
| block diversity cost | Float   | b          | beta             | energy cost for interblock diversity due to alignment merger                                                  |
| circular genomes     | Boolean | c          | circular         | toggle if input genomes are circular                                                                          |
| pairwise sensitivity | String  | s          | sensitivity      | controls the pairwise genome alignment sensitivity of minimap 2. Currently only accepts "5", "10" or "20"     |
| maximum self-maps    | Integer | x          | max-self-map     | maximum number of iterations to perform block self maps per pairwise graph merger                             |
| enforce uppercase    | Boolean | u          | upper-case       | toggle to force genomes to uppercase characters                                                               |
| distance calculator  | String  | d          | distance-backend | only accepts "native" or "mash"                                                                               |
| alignment kernel     | String  | k          | alignment-kernel | only accepts "minimap2" or "mmseqs"                                                                           |
| kmer length (mmseqs) | Integer | K          | kmer-length      | kmer length, only used for mmseqs2 alignment kernel. If not specified will use mmseqs default.                |
| consistency check    | Boolean | t          | test             | toggle to activate consistency check: verifies that input genomes can be exactly reconstructed from the graph |

## Arguments
Expects one or more fasta files.
Multiple records within one file are treated as separate genomes
Fasta files can be optionally gzipped.

## Output
Prints the constructed pangraph as a JSON to _stdout_.

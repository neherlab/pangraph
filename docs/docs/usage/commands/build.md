# Command: build

## Introduction

Vivamus luctus egestas leo. Etiam urna velit, aliquam in vulputate at, scelerisque vitae mi. Pellentesque ac magna purus. Nunc fermentum tortor ac porta dapibus. In rutrum ac purus sit amet tempus. Aliquam erat volutpat. Sed elementum, eros sed suscipit adipiscing, dolor diam aliquet lectus, vitae dictum justo urna id quam.

## Context

Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo. Nemo enim ipsam voluptatem quia voluptas sit aspernatur aut odit aut fugit, sed quia consequuntur magni dolores eos qui ratione voluptatem sequi nesciunt.

## Arguments

Neque porro quisquam est, qui dolorem ipsum quia dolor sit amet, consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt ut labore et dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi consequatur.

| Argument     | Description                                                                                                                    | Additional Information                                                                                                                                                                                                                                         |
|--------------|--------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| INPUT_FASTAS | Path(s) to zero, one, or multiple FASTA files with input sequences. Each record within a file is treated as a separate genome. | - Supports plain or compressed FASTA files.<br/> - Automatic decompression for `gz`, `bz2`, `xz`, `zstd` formats, based on file extension.<br/> - Reads from stdin if no files are provided.<br/> - [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) |

## Options

| Option                            | Description                                           | Default | Details                                                                                         |
|-----------------------------------|-------------------------------------------------------|---------|-------------------------------------------------------------------------------------------------|
| -o, --output-json \<OUTPUT_JSON\> | Path to output JSON file with the resulting pangraph. | `-`     | - Supports output compression.<br/> - Creates directories as needed.<br/> - Use "-" for stdout. |
| -h, --help                        | Print help information.                               |         | Use '-h' for a summary.                                                                         |

## Alignment

| Option                                      | Description                                                                          | Default        | Details                                                   |
|---------------------------------------------|--------------------------------------------------------------------------------------|----------------|-----------------------------------------------------------|
| -l, --len \<INDEL_LEN_THRESHOLD\>           | Minimum block size for the alignment graph (in nucleotides).                         | `100`          |                                                           |
| -a, --alpha \<ALPHA\>                       | Energy cost for introducing a junction due to alignment merger.                      | `100`          |                                                           |
| -b, --beta \<BETA\>                         | Energy cost for interblock diversity due to alignment merger.                        | `10`           |                                                           |
| -s, --sensitivity \<SENSITIVITY\>           | Sets pairwise alignment sensitivity.                                                 | `10`           |                                                           |
| -K, --kmer-length \<KMER_LENGTH\>           | Sets kmer length for the MMseqs2 aligner.                                            |                |                                                           |
| -c, --circular                              | Toggle if input genomes are circular.                                                |                |                                                           |
| -u, --upper-case                            | Transforms all sequences to upper case.                                              |                |                                                           |
| -x, --max-self-map \<MAX_SELF_MAP\>         | Maximum number of alignment rounds per pairwise graph merger.                        | `100`          |                                                           |
| -d, --distance-backend \<DISTANCE_BACKEND\> | Backend for genome similarity estimation.                                            | `native`       | Possible values: `native`, `mash`                         |
| -k, --alignment-kernel \<ALIGNMENT_KERNEL\> | Backend for pairwise genome alignment.                                               | `minimap2-lib` | Possible values: `minimap2-lib`, `minimap2-cli`, `mmseqs` |
| -f, --verify                                | Verify that the original sequences can be reconstructed from the resulting pangraph. |                |                                                           |
| --seed \<SEED\>                             | Random seed for block ID generation.                                                 |                |                                                           |

## Verbosity

| Option                    | Description                            | Default | Details                                                           |
|---------------------------|----------------------------------------|---------|-------------------------------------------------------------------|
| -j, --jobs \<JOBS\>       | Number of processing jobs.             |         | Uses all CPU threads if not specified.                            |
| --verbosity \<VERBOSITY\> | Set verbosity level of console output. | `warn`  | Possible values: `off`, `error`, `warn`, `info`, `debug`, `trace` |
| --silent                  | Disable all console output.            |         | Equivalent to `--verbosity=off`.                                  |
| -v, --verbose...          | Increase verbosity of console output.  |         | Add occurrences for higher verbosity.                             |
| -q, --quiet...            | Reduce verbosity of console output.    |         | Add occurrences for quieter output.                               |

### Examples

#### Help

Display command line options for `build` command

```bash
pangraph build --help

```

#### Basic (to stdout)

Take `my.fasta` and print Pangraph JSON to standard output

```bash
pangraph build my.fasta
```

#### Basic (to file)

Take `my.fasta` and output Pangraph JSON to a file

```bash
pangraph build my.fasta -o my.pangraph.json
```


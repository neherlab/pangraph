# Export

## Description
Export a pangraph to a chosen file format(s)

## Options
| Name                | Type    | Short Flag | Long Flag           | Description                                                                       |
| :------------------ | :------ | :--------- | :------------------ | :-------------------------------------------------------------------------------- |
| Edge minimum length | Integer | ell        | edge-minimum-length | blocks below this length cutoff will be ignored for edges in graph                |
| Edge maximum length | Integer | elu        | edge-maximum-length | blocks above this length cutoff will be ignored for edges in graph                |
| Edge minimum depth  | Integer | edl        | edge-minimum-depth  | blocks below this depth cutoff will be ignored for edges in graph                 |
| Edge maximum depth  | Integer | edu        | edge-maximum-depth  | blocks above this depth cutoff will be ignored for edges in graph                 |
| Minimum length      | Integer | ll         | minimum-length      | blocks below this length cutoff will be ignored for export                        |
| Maximum length      | Integer | lu         | maximum-length      | blocks above this length cutoff will be ignored for export                        |
| Minimum depth       | Integer | dl         | minimum-depth       | blocks below this depth cutoff will be ignored for export                         |
| Maximum depth       | Integer | du         | maximum-depth       | blocks above this depth cutoff will be ignored for export                         |
| No duplications     | Boolean | nd         | no-duplications     | do not export any block that contains at least one strain repeated more than once |
| Output directory    | String  | o          | output-directory    | path to directory where output will be stored (default: `export`)                 |
| Prefix              | String  | p          | prefix              | basename of exported files (default: `pangraph`)                                  |
| GFA                 | Boolean | ng         | no-export-gfa       | toggles whether pangraph is exported as GFA.                                      |
| PanX                | Boolean | pX         | export-panX         | toggles whether pangraph is exported to panX visualization compatible format.     |

## Arguments
Zero or one pangraph file which must be formatted as a JSON.
If no file path is given, reads from _stdin_.
In either case, the stream can be optionally gzipped.

## Output
Outputs the constructed pangraph to the selected formats at the user-supplied paths.

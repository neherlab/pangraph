# Polish

## Description
Realigns blocks of a multiple sequence alignment pangraph with an external multiple sequence alignment tool.

## Options
Name | Type | Short Flag | Long Flag | Description
:-------------- | :------- | :------ | :------- | :-------------------------
Maximum Length | Integer | l | length | cutoff above which the block is not realigned externally
Preserve Case | Bool | c | preserve-case | ensure case (upper/lower) is preserved after realignment

## Arguments
Zero or one pangraph file which must be formatted as a JSON.
If no file path is given, reads from _stdin_.
In either case, the stream can be optionally gzipped.

## Output
Outputs the polished pangraph to _stdout_.

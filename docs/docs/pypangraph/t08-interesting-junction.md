---
sidebar_position: 10
---

# A look at an interesting junction

Once the summary statistics from the [previous part](t07-junction-stats.md) have pointed out a junction worth investigating, the natural next step is to **zoom in on that one junction** and look at it in detail.

`BackboneJunctions` offers two methods for this:

- **`positions()`** — where each occurrence of the junction sits on the genomes that carry it. Useful e.g. for cross-referencing with annotation files.
- **`sequences()`** — the actual DNA spanning the junction on every genome, returned as `Bio.SeqRecord` objects ready to be written to FASTA for further analysis (multiple sequence alignment, secondary pangraph construction, BLAST searches, ...).


## The example: a candidate IS insertion

As a running example we will use the same junction introduced [in part 6](t06-junctions-intro.md): the edge `3156970751805415521_f__4335229004353524956_f`. Looking up its row in the statistics dataframe:

```python
import pypangraph as pp

graph = pp.Pangraph.from_json("staph.json.gz")
bj = pp.junctions.BackboneJunctions(graph, L_thr=500)
stats = bj.stats()

edge = "3156970751805415521_f__4335229004353524956_f"
print(stats.loc[edge, ["n_isolates", "n_non_empty", "n_categories", "accessory_length"]])
# n_isolates            15
# n_non_empty            3
# n_categories           2
# accessory_length    1512
```

The junction is universal (`n_isolates == 15`), only **3 of the 15 isolates** carry accessory content between the flanking core blocks (`n_non_empty == 3`), and there are only two distinct accessory variants: empty, and a single shared insert of about **1.5 kb** (`accessory_length == 1512`). This is the textbook profile of a recently acquired **insertion sequence** present in a minority of the population.

## Locating the junction on each genome

The `positions()` method returns a `pandas.DataFrame` indexed by `(isolate, edge)`, with the genomic coordinates of the two flanking core blocks plus a strand flag. We slice on our edge of interest with `.xs`:

```python
pos = bj.positions().xs(edge, level=1)
print(pos.head(6))
#                left_start  left_end  right_start  right_end  strand
# iso
# NZ_CP132362.1       75741     83703        83703      89270   False
# NZ_LR822061.1     1725093   1733055      1733055    1738656   False
# NZ_CP077852.1     1234589   1240190      1240190    1248152    True
# NZ_CP162433.1     1568982   1576944      1578456    1584065   False
# NZ_CP034011.1      414663    420264       420264     428226    True
# NZ_CP092558.1      392241    397850       397850     405804    True
```

The columns are:

- **`left_start`**, **`left_end`** — genomic coordinates of the left flanking core block.
- **`right_start`**, **`right_end`** — genomic coordinates of the right flanking core block.
- **`strand`** — `True` if the junction appears in canonical edge orientation on this genome, `False` if reverse-complemented.

The accessory content of the junction sits in the interval `[left_end, right_start)`. Notice how on most isolates `left_end == right_start` (the two flanks sit directly adjacent: empty junction), but on `NZ_CP162433.1` the gap is `1578456 - 1576944 = 1512` bp — the IS insert. Hunting for rows with `right_start - left_end > 0` is a quick way to find the IS carriers:

```python
pos.eval("right_start - left_end").rename("insert_len").to_frame().query("insert_len > 0")
#                insert_len
# iso
# NZ_CP162433.1        1512
# NZ_AP024511.1        1520
# NZ_CP080550.1        1515
```

The small per-isolate differences (1512 / 1515 / 1520 bp) come from indels in the IS sequence; the **consensus** length of the shared accessory block — which is what `accessory_length` in the stats reports — is exactly 1512 bp.

:::info circular genomes

For circular genomes the coordinates may wrap around the origin, in which case `left_start > left_end` on the isolate where the junction straddles position 0.

:::

<details>
    <summary>cross-referencing with genome annotations</summary>

    Knowing where each junction sits on each genome lets us cross-reference these coordinates with annotation files (e.g. GFF) — for instance to check whether the IS insertion overlaps any predicted CDS, prophage, or transposase annotation on a given isolate.

    The sketch below picks one IS-carrying isolate and asks for annotations overlapping the insert interval:

    ```python
    iso = "NZ_CP162433.1"
    insert_start = pos.loc[iso, "left_end"]
    insert_end = pos.loc[iso, "right_start"]

    # `annotations` is a hypothetical dataframe with columns
    # "iso", "start", "end", "feature" loaded from a GFF file
    overlaps = annotations.query(
        "iso == @iso and start < @insert_end and end > @insert_start"
    )
    print(overlaps[["start", "end", "feature"]])
    ```

    The same logic can be applied to every IS-carrying isolate in a loop to confirm that the insert is annotated as a mobile element on each of them.

</details>

## Extracting the junction sequences

The `sequences()` method returns one `Bio.SeqRecord` per isolate, spanning **left flank + accessory center + right flank**, all co-oriented to the canonical edge direction:

```python
records = bj.sequences(edge)
for r in records:
    print(f"  {r.id}: {len(r.seq)} bp")
#   NZ_CP132362.1: 13529 bp
#   NZ_LR822061.1: 13563 bp
#   NZ_CP077852.1: 13563 bp
#   NZ_CP162433.1: 15083 bp  <-- IS carrier
#   NZ_CP034011.1: 13563 bp
#   NZ_CP092558.1: 13563 bp
#   NZ_CP062358.1: 13563 bp
#   NZ_CP132372.1: 13563 bp
#   NZ_AP024511.1: 15083 bp  <-- IS carrier
#   NZ_CP181145.1: 13563 bp
#   NZ_CP157420.1: 13563 bp
#   NZ_CP080550.1: 15083 bp  <-- IS carrier
#   NZ_CP169947.1: 13563 bp
#   NZ_CP022905.1: 13563 bp
#   NZ_CP035791.1: 13563 bp
```

The three IS carriers stand out as ~1.5 kb longer than the rest. Each record has:

- **`id`** — the isolate name.
- **`description`** — the canonical edge string.
- **`seq`** — the DNA from the start of the left flank to the end of the right flank.

<details>
    <summary>co-orientation guarantees</summary>

    All returned sequences are co-oriented to the canonical edge: sequences from genomes where the junction is found in the inverted orientation are automatically reverse-complemented before being returned. Individual accessory blocks within each sequence are also placed in the correct orientation based on their strand. This makes the flanking core blocks directly alignable across isolates without further pre-processing.

</details>

### Saving to a FASTA file

The records are ready to be written out with `Bio.SeqIO`:

```python
from Bio import SeqIO

SeqIO.write(records, f"junction_{edge}.fasta", "fasta")
```

From there the FASTA can be fed to any downstream tool — for example a multiple sequence alignment with [MAFFT](https://mafft.cbrc.jp/alignment/software/):

```bash
mafft junction_3156970751805415521_f__4335229004353524956_f.fasta > junction_aligned.fasta
```

Or a second, junction-specific pangraph build to resolve sub-block structure in the IS region:

```bash
pangraph build -l 100 -s 20 junction_aligned.fasta -o junction.json.gz
```

### Trimming away the flanks

If only the accessory portion of each junction is of interest, the flank lengths from the statistics dataframe can be used to slice off the core flanks:

```python
left_len = stats.loc[edge, "left_core_length"]   # 5601
right_len = stats.loc[edge, "right_core_length"] # 7962

for r in records:
    accessory = r.seq[left_len:-right_len] if right_len > 0 else r.seq[left_len:]
    print(f"  {r.id}: accessory = {len(accessory)} bp")
#   NZ_CP132362.1: accessory =    0 bp
#   NZ_LR822061.1: accessory =    0 bp
#   ...
#   NZ_CP162433.1: accessory = 1520 bp  <-- the IS
#   ...
```

Empty junctions cleanly trim to `0` bp; the IS carriers each retain ~1.5 kb of accessory content as expected.

## Exporting many junctions at once

<details>
    <summary>batch export of all non-transitive junctions to FASTA files</summary>

    The same workflow scales naturally to every non-transitive junction in the graph: just iterate over the relevant rows of the statistics dataframe and call `sequences()` on each edge.

    ```python
    import os

    os.makedirs("junctions", exist_ok=True)

    non_transitive = stats[~stats["is_transitive"]].index
    for e in non_transitive:
        recs = bj.sequences(e)
        if recs:
            SeqIO.write(recs, f"junctions/{e}.fasta", "fasta")

    print(f"Exported {len(non_transitive)} junction files")
    # Exported 132 junction files
    ```

</details>

## An exercise for the reader

#TODO: give another example of an interesting junction to look at
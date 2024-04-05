# class specifications

## Hit

- name: str - sequence id
- length: int - sequence length
- start: int - mapping start position (0-based indexing)
- stop: int - mapping end position (0-based indexing, right-end excluded, same as minimap paf output)
- seq: str - optional, sequence. **removed**.

## Alignment

- qry: Hit - hit on the qry sequence
- reff: Hit - hit on the reference sequence
- matches: int - number of matches on the cigar string
- length: int - total length of the alignment
- quality: int - match quality
- orientation: str - fwd or rev match
- cigar: str - cigar string
- divergence: float - approximate sequence divergence
- align: float

# notes

## minimap 2 paf files

- start/end coordinates are given in the same convention of Python array slicing:
  - 0-based indexing
  - right-end excluded

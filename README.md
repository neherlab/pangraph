# pypangraph - rust version

upgrade of the `pypangraph` library to parse

## dev roadmap

Install the package in dev mode:

```bash
pip install -e .
```

- [x] parse the graph
  - [x] update alignment format
- [x] produce block statistics
- [ ] produce alignments
  - [ ] set of unaligned fasta records
  - [ ] biopython alignment, but without insertions
  - [ ] biopython alignment, without insertions or deletions
- [ ] core-sinteny analysis
- [ ] find core-junctions
- [ ] run core-junction analysis
  - [ ] rebuild sub-graphs
  - [ ] analyze sub-graphs
- [ ] locator: from position on the genome to position in the graph, and vice-versa.
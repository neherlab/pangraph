# The structure of Pangraph output file

- aim of the tutorial: more details on json output file
    - path object
    - block object
- frequency distribution of blocks? Weighted by length?
- 

## The structure of `pangraph.json`


### block entries

```json
{
    "id": "KZJIDOXBAV",
    "sequence": "AAGGTGGGTAATCATTTTGATAAGTGAT...",
    "gaps": {...},
    "mutate": [...],
    "insert": [...],
    "delete": [...],
    "positions": [...]
}
```

### path entries

```json
{
    "name": "NZ_CP010242",
    "blocks": [
        { "id": "NFTNKNMFIC", ... },
        { "id": "YTLSRRNNGL", ... },
        { "id": "HDOKGYMDQR", ... }
    ]
    ...
}
```


"""
Example demonstrating how to read Pangraph JSON using Python classes generated from JSON schema.
See README.md in the parent directory for instructions.
```

"""
import sys
import json
from dacite import from_dict
from Pangraph_model import Pangraph


def read_pangraph(filepath: str = None) -> Pangraph:
    """
    Read pangraph JSON and convert it to Python dataclasses (recursively). If no filename is provided, reads from stdin.
    """
    source = sys.stdin if filepath is None else open(filepath)
    with source as f:
        json_data = json.load(f)
        return from_dict(data_class=Pangraph, data=json_data)


if __name__ == "__main__":
    filepath = sys.argv[1] if len(sys.argv) > 1 else None

    # pangraph = parse_file_as(Pangraph, filepath)
    pangraph = read_pangraph(filepath)

    for bid, block in list(pangraph.blocks.items())[0:2]:
        print(f"{bid=}")
        print(f"{block.consensus=}")
        print(f"{list(block.alignments.values())[0]=}")

    print("")

    for nid, node in list(pangraph.nodes.items())[0:3]:
        print(f"{nid=}")
        print(f"{node=}")

    print("")

    for pid, paths in list(pangraph.paths.items())[0:3]:
        print(f"{pid=}")
        print(f"{paths=}")

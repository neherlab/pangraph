import json
import os
import sys
import importlib.resources as pkg_resources
from os.path import join, dirname
from jsonschema import validate
from typing import Any
from Pangraph_model import Pangraph
from dacite import from_dict


def read_pangraph(filepath: str = None) -> Pangraph:
    """
    Read pangraph JSON, validate it, and convert it to Python dataclasses (recursively).
    If no filename is provided, reads from stdin.
    """
    source = sys.stdin if filepath is None else open(filepath)
    with source as f:
        json_data = json.load(f)
        validate_pangraph({"pangraph": json_data}) # Extra field for the dummy container
        return from_dict(data_class=Pangraph, data=json_data)


def validate_pangraph(data: Any):
    """
    Validate pangraph json data against JSON schema
    """
    schema = load_schema()
    validate(instance=data, schema=schema)


def load_schema():
    """
    Load json schema from the package resources
    """
    if is_dev_mode():
        with open(join(dirname(os.path.abspath(__file__)), 'Pangraph.schema.json'), "r") as f:
            return json.load(f)
    else:
        with pkg_resources.open_text(__name__, 'Pangraph.schema.json') as f:
            return json.load(f)


def is_dev_mode():
    return 'site-packages' not in os.path.abspath(__file__)


# Remove this
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

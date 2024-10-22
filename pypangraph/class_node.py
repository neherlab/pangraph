def parse_strandedness(strand: str) -> bool:
    """Parse strandedness from a string to a boolean:
    + -> True,
    - -> False
    """
    match strand:
        case "+":
            return True
        case "-":
            return False
        case _:
            raise ValueError(f"Strand {strand} not recognized.")


class Node:
    """Pangraph node object. It has attributes:
    - name (str): strain name
    - circular (bool): whether the path object is circular
    - nodes (list): list of node ids in the path
    - nuc_len (int): total length of the path in base-pairs
    """

    def __init__(self, node):
        self.id = node["id"]
        self.block_id = node["block_id"]
        self.path_id = node["path_id"]
        self.strand = parse_strandedness(node["strand"])
        self.position = tuple(node["position"])

    def __str__(self):
        return f"node {self.id} (block={self.block_id}, path={self.path_id})"

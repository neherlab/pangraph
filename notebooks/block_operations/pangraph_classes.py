class Insertion:
    def __init__(self, pos: int, ins: str):
        self.pos = pos
        self.ins = ins

    def __str__(self):
        return f"Insertion(pos={self.pos}, ins={self.ins})"


class Deletion:
    def __init__(self, pos: int, length: int):
        self.pos = pos
        self.length = length

    def __str__(self):
        return f"Deletion(pos={self.pos}, length={self.length})"


class Substitution:
    def __init__(self, pos: int, alt: str):
        self.pos = pos
        self.alt = alt

    def __str__(self):
        return f"Substitution(pos={self.pos}, alt={self.alt})"


class Edits:
    """
    any query sequence can be transformed into the reference sequence by applying
    a set of edits.
    """

    def __init__(
        self,
        ins: list[Insertion] = [],
        dels: list[Deletion] = [],
        subs: list[Substitution] = [],
    ):
        self.ins = ins
        self.dels = dels
        self.subs = subs

    def __str__(self):
        msg = "Edits("
        for i in self.ins:
            msg += f"\nIns -> {i}\t"
        for d in self.dels:
            msg += f"\nDel -> {d}\t"
        for s in self.subs:
            msg += f"\nSub -> {s}\t"
        msg += "\n)"
        return msg


class Node:
    def __init__(
        self,
        node_id: str,
        block_id: str,
        path_id: str,
        position: tuple[int, int],
        strand: bool,
    ):
        self.node_id = node_id
        self.block_id = block_id
        self.path_id = path_id
        self.position = position
        self.strand = strand

    def __str__(self):
        return f"Node(node_id={self.node_id}, block_id={self.block_id}, path_id={self.path_id}, position={self.position}, strand={self.strand})"


class Block:
    def __init__(self, block_id: str, sequence: str, alignment: dict[str, Edits]):
        self.id = block_id
        self.consensus = sequence
        self.alignment = alignment

    def __str__(self):
        msg = f"Block(\nblock_id={self.id}\nconsensus={self.consensus}\nalignment:"
        for k, v in self.alignment.items():
            msg += f"\n{k} -> {v}"
        msg += "\n)"
        return msg


class Path:
    def __init__(
        self,
        path_id: str,
        nodes: list[Node],
        tot_len: int,
        offset: int,
        circular: bool,
        positions: dict[str, int],
    ):
        self.path_id = path_id
        self.nodes = nodes
        self.tot_len = tot_len
        self.offset = offset
        self.circular = circular
        self.positions = positions

    def __str__(self):
        return f"Path(path_id={self.path_id}, nodes={self.nodes}, tot_len={self.tot_len}, offset={self.offset}, circular={self.circular}, positions={self.positions})"


class Pangraph:
    def __init__(self, blocks: dict[str, Block], paths: dict[str, Path]):
        self.blocks = blocks
        self.paths = paths

    def __str__(self):
        return f"PanGraph(blocks={self.blocks}, paths={self.paths})"

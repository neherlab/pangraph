from dataclasses import dataclass, field
from typing import List, Tuple, Dict


@dataclass
class SimpleNode:
    bid: int  # block id
    strand: bool  # strand

    def __invert__(self) -> "SimpleNode":
        return SimpleNode(self.bid, not self.strand)

    def __eq__(self, other: object) -> bool:
        return self.bid == other.bid and self.strand == other.strand

    def __hash__(self) -> int:
        return hash((self.bid, self.strand))

    def __repr__(self) -> str:
        s = "+" if self.strand else "-"
        return f"({self.bid}{s})"


@dataclass
class Edge:
    n1: SimpleNode
    n2: SimpleNode

    def __invert__(self) -> "Edge":
        return Edge(~self.n2, ~self.n1)

    def oriented_equal(self, other: "Edge") -> bool:
        return self.n1 == other.n1 and self.n2 == other.n2

    def __eq__(self, other: object) -> bool:
        return self.oriented_equal(other) or self.oriented_equal(~other)

    def side_hash(self) -> int:
        return hash((self.n1, self.n2))

    def __hash__(self) -> int:
        return self.side_hash() ^ (~self).side_hash()

    def __repr__(self) -> str:
        return f"[{self.n1}|{self.n2}]"


@dataclass
class SimplePath:
    nodes: List[SimpleNode]
    circular: bool

    def to_edges(self) -> List[Edge]:
        edges = [
            Edge(self.nodes[i], self.nodes[i + 1]) for i in range(len(self.nodes) - 1)
        ]
        if self.circular:
            edges.append(Edge(self.nodes[-1], self.nodes[0]))
        return edges


@dataclass
class GFA:
    """
    Class to represent a GFA object.
    """

    header: str = "H\tVN:Z:1.0"
    segments: List[Tuple[str, str, int, int, bool]] = field(
        default_factory=list
    )  # List of (name, sequence, read_count, length, duplicated)
    links: Dict[Edge, int] = field(default_factory=dict)  # Map of Edge to read count
    paths: List[Tuple[str, List[str], List[str], bool]] = field(
        default_factory=list
    )  # List of (path_name, segments, overlaps, circular)

    def add_segment(
        self,
        name: str,
        sequence: str = "*",
        read_count: int = 0,
        length: int = 0,
        duplicated: bool = False,
    ) -> None:
        self.segments.append((name, sequence, read_count, length, duplicated))

    def add_link(self, edge: Edge) -> None:
        if edge in self.links:
            self.links[edge] += 1
        else:
            self.links[edge] = 1

    def add_path(
        self,
        path_name: str,
        segments: List[str],
        overlaps: List[str],
        circular: bool = False,
    ) -> None:
        self.paths.append((path_name, segments, overlaps, circular))


def convert_pangraph_to_simple_paths(
    pangraph, minimal_length: int, export_duplicated: bool
) -> List[SimplePath]:
    """
    Converts a Pangraph object to a list of SimplePath objects.

    Args:
        pangraph (Pangraph): The Pangraph object.
        minimal_length (int): Minimum length of block consensus sequence to include.
        export_duplicated (bool): Whether to include duplicated blocks.

    Returns:
        List[SimplePath]: A list of SimplePath objects.
    """
    simple_paths = []

    for path in pangraph.paths.values():
        nodes = []
        for node_id in path.nodes:
            node = pangraph.nodes[node_id]
            block = pangraph.blocks[node.block_id]

            # Filter by minimal length
            if len(block.consensus) < minimal_length:
                continue

            # Check duplication if needed
            if not export_duplicated:
                is_duplicated = any(
                    path.nodes.count(nid) > 1 for nid in block.alignments.keys()
                )
                if is_duplicated:
                    continue

            nodes.append(SimpleNode(node.block_id, node.strand))

        # Add path if nodes exist
        if nodes:
            simple_paths.append(SimplePath(nodes=nodes, circular=path.circular))

    return simple_paths


def convert_simple_paths_to_gfa(
    simple_paths: List[SimplePath], block_data: Dict[int, Tuple[str, int]]
) -> GFA:
    """
    Converts a list of SimplePath objects and block data to a GFA object.

    Args:
        simple_paths (List[SimplePath]): The list of SimplePath objects.
        block_data (Dict[int, Tuple[str, int]]): Block data as {block_id: (sequence, length)}.

    Returns:
        GFA: The GFA object.
    """
    gfa = GFA()

    # Add segments
    for block_id, (sequence, length) in block_data.items():
        segment_name = f"block_{block_id}"
        duplicated = any(
            path.nodes.count(SimpleNode(block_id, True))
            + path.nodes.count(SimpleNode(block_id, False))
            > 1
            for path in simple_paths
        )
        gfa.add_segment(segment_name, sequence, 0, length, duplicated)

    # Add links and paths
    for i, path in enumerate(simple_paths):
        edges = path.to_edges()
        for edge in edges:
            gfa.add_link(edge)
        segments = [
            f"block_{node.bid}{'+' if node.strand else '-'}" for node in path.nodes
        ]
        overlaps = ["*"] * (len(segments) - 1)
        gfa.add_path(f"path_{i}", segments, overlaps, path.circular)

    return gfa


def write_gfa_to_file(gfa: GFA, output_file: str) -> None:
    """
    Writes a GFA object to a file.

    Args:
        gfa (GFA): The GFA object.
        output_file (str): Path to the output GFA file.
    """
    with open(output_file, "w") as f:
        # Write header
        f.write(f"{gfa.header}\n")

        # Write segments
        f.write("# blocks\n")
        for name, sequence, read_count, length, duplicated in gfa.segments:
            duplicated_tag = "\tDP:Z:duplicated" if duplicated else ""
            f.write(
                f"S\t{name}\t{sequence}\tRC:i:{read_count}\tLN:i:{length}{duplicated_tag}\n"
            )

        # Write links
        f.write("# edges\n")
        for edge, read_count in gfa.links.items():
            n1_orientation = "+" if edge.n1.strand else "-"
            n2_orientation = "+" if edge.n2.strand else "-"
            f.write(
                f"L\tblock_{edge.n1.bid}\t{n1_orientation}\tblock_{edge.n2.bid}\t{n2_orientation}\t*\tRC:i:{read_count}\n"
            )

        # Write paths
        f.write("# paths\n")
        for path_name, segments, overlaps, circular in gfa.paths:
            circular_tag = "\tTP:Z:circular" if circular else ""
            f.write(
                f"P\t{path_name}\t{','.join(segments)}\t{','.join(overlaps)}{circular_tag}\n"
            )


# Example usage:
# simple_paths = [SimplePath(...), ...]
# block_data = {block_id: (sequence, length), ...}
# gfa = convert_simple_paths_to_gfa(simple_paths, block_data)
# write_gfa_to_file(gfa, "output.gfa")

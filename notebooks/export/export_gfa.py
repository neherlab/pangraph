from dataclasses import dataclass, field
from typing import List, Tuple, Dict
from .circularize_utils import SimpleNode, Edge
from .utils import Pangraph
from collections import defaultdict


@dataclass
class GFA_segment:
    name: str  # block id
    sequence: str  # block sequence
    depth: int  # number of block occurrences
    length: int  # block length
    duplicated: bool  # whether the block is duplicated


@dataclass
class GFA_links:
    edge_ct: defaultdict[Edge, int] = field(default_factory=lambda: defaultdict(int))

    def add_edge(self, bid1: int, strand1: bool, bid2: int, strand2: bool):
        n1 = SimpleNode(bid1, strand1)
        n2 = SimpleNode(bid2, strand2)
        edge = Edge(n1, n2)
        self.edge_ct[edge] += 1


@dataclass
class GFA_path:
    path_name: str
    # List of block ids and strandedness
    segments: List[SimpleNode]
    circular: bool


@dataclass
class GFA:
    """
    Class to represent a GFA object.
    """

    header: str = "H\tVN:Z:1.0"
    # Dictionary id -> segment
    segments: Dict[int, GFA_segment] = field(default_factory=dict)
    # Dictionary edge -> read count
    links: GFA_links = field(default_factory=GFA_links)
    # List of paths
    paths: List[GFA_path] = field(default_factory=list)

    def to_file(self, fname: str, export_sequence: bool) -> None:
        """
        Write GFA object to file.
        """
        with open(fname, "w") as f:
            # Write header
            f.write(f"{self.header}\n")

            # Write segments
            f.write("# blocks\n")
            for segment in self.segments.values():
                segment_seq = f"{segment.sequence}" if export_sequence else "*"
                duplicated_tag = "\tDP:Z:duplicated" if segment.duplicated else ""
                read_count = segment.depth * segment.length
                f.write(
                    f"S\t{segment.name}\t{segment_seq}\tRC:i:{read_count}\tLN:i:{segment.length}{duplicated_tag}\n"
                )

            # Write links
            f.write("# edges\n")
            for edge, read_count in self.links.edge_ct.items():
                n1_orientation = "+" if edge.n1.strand else "-"
                n2_orientation = "+" if edge.n2.strand else "-"
                f.write(
                    f"L\t{edge.n1.bid}\t{n1_orientation}\t{edge.n2.bid}\t{n2_orientation}\t*\tRC:i:{read_count}\n"
                )

            # Write paths
            f.write("# paths\n")
            for path in self.paths:
                circular_tag = "\tTP:Z:circular" if path.circular else ""
                path_segments = [
                    f"{node.bid}{('+' if node.strand else '-')}"
                    for node in path.segments
                ]
                f.write(
                    f"P\t{path.path_name}\t{','.join(path_segments)}\t*{circular_tag}\n"
                )


def graph_to_gfa_segment_dict(graph: Pangraph) -> Dict[int, GFA_segment]:
    """
    Convert the set of blocks of a pangraph to a dictionary of GFA segments,
    whose keys are the block ids.
    """
    segments = {}
    for block_id, block in graph.blocks.items():
        isolates = [graph.nodes[node_id].path_id for node_id in block.alignments.keys()]
        is_duplicated = len(set(isolates)) < len(isolates)
        segments[block_id] = GFA_segment(
            name=block_id,
            sequence=block.consensus,
            depth=len(block.alignments),
            length=len(block.consensus),
            duplicated=is_duplicated,
        )
    return segments


def graph_to_gfa_paths(graph: Pangraph) -> List[GFA_path]:
    """
    Convert the paths of a pangraph to a list of GFA paths.
    """
    paths = []
    for path_id, path in graph.paths.items():
        segments = []
        for node_id in path.nodes:
            node = graph.nodes[node_id]
            segments.append(SimpleNode(bid=node.block_id, strand=node.strand))
        paths.append(
            GFA_path(path_name=path.name, segments=segments, circular=path.circular)
        )
    return paths


def filter_gfa_path(
    path: GFA_path,
    segments: Dict[int, GFA_segment],
    min_length: int,
    min_depth: int,
    export_duplicated: bool,
) -> List[GFA_path]:
    """
    Optionally filters out segments in a path that are:
    - shorter than min_length
    - present less than min_depth times
    - duplicated
    """
    new_segments = []
    for node in path.segments:
        if min_length is not None and segments[node.bid].length < min_length:
            continue
        if min_depth is not None and segments[node.bid].depth < min_depth:
            continue
        if not export_duplicated and segments[node.bid].duplicated:
            continue
        new_segments.append(node)
    return GFA_path(
        path_name=path.path_name, segments=new_segments, circular=path.circular
    )


def gfa_links_from_paths(paths: List[GFA_path]) -> GFA_links:
    links = GFA_links()
    for path in paths:
        for i in range(len(path.segments) - 1):
            seg1 = path.segments[i]
            seg2 = path.segments[i + 1]
            links.add_edge(seg1.bid, seg1.strand, seg2.bid, seg2.strand)
        if path.circular:
            seg1 = path.segments[-1]
            seg2 = path.segments[0]
            links.add_edge(seg1.bid, seg1.strand, seg2.bid, seg2.strand)
    return links


def graph_to_gfa(
    graph: Pangraph,
    min_len: int = None,
    min_depth: int = None,
    export_duplicated: bool = True,
) -> GFA:
    """
    Convert a pangraph to a GFA object.
    Filter out blocks that are shorter than min_len or are duplicated.
    """
    GFA_segments = graph_to_gfa_segment_dict(graph)
    GFA_paths = graph_to_gfa_paths(graph)
    GFA_paths = [
        filter_gfa_path(
            path,
            GFA_segments,
            min_depth=min_depth,
            min_length=min_len,
            export_duplicated=export_duplicated,
        )
        for path in GFA_paths
    ]
    # remove empty paths
    GFA_paths = [path for path in GFA_paths if len(path.segments) > 0]

    # filter GFA segments that are not in any path
    segment_ids = set()
    for path in GFA_paths:
        for segment in path.segments:
            segment_ids.add(segment.bid)
    GFA_segments = {k: v for k, v in GFA_segments.items() if k in segment_ids}

    GFA_links = gfa_links_from_paths(GFA_paths)

    return GFA(segments=GFA_segments, links=GFA_links, paths=GFA_paths)


def export_graph_to_gfa(
    graph: Pangraph,
    output_file: str,
    min_len: int,
    min_depth: int,
    export_duplicated: bool,
    export_sequence: bool = False,
) -> None:
    """
    Export a pangraph to a GFA file.
    """
    gfa = graph_to_gfa(
        graph,
        min_len=min_len,
        min_depth=min_depth,
        export_duplicated=export_duplicated,
    )

    gfa.to_file(output_file, export_sequence)

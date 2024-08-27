from collections import defaultdict, Counter


class Node:
    """Combination of block id and strandedness"""

    def __init__(self, bid: str, strand: bool) -> None:
        self.id = bid
        self.strand = strand

    def invert(self) -> "Node":
        return Node(self.id, not self.strand)

    def __invert__(self) -> "Node":
        return self.invert()

    def __eq__(self, other: object) -> bool:
        return self.id == other.id and self.strand == other.strand

    def __hash__(self) -> int:
        return hash((self.id, self.strand))

    def __repr__(self) -> str:
        s = "+" if self.strand else "-"
        return f"[{self.id}|{s}]"

    def to_str_id(self):
        s = "f" if self.strand else "r"
        return f"{self.id}_{s}"

    @staticmethod
    def from_str_id(t) -> "Node":
        bid = t.split("_")[0]
        strand = True if t.split("_")[1] == "f" else False
        return Node(bid, strand)


class Path:
    """A path is a list of nodes"""

    def __init__(self, nodes=[], circular=None) -> None:
        self.nodes = nodes
        self.circular = circular

    def add_left(self, node: Node) -> None:
        self.nodes.insert(0, node)

    def add_right(self, node: Node) -> None:
        self.nodes.append(node)

    def rotate_to(self, bid: str, strand: bool) -> "Path":
        assert self.circular, "Path is not circular"
        assert bid in [n.id for n in self.nodes], "Block not in path"
        n = Node(bid, strand)
        if n in self.nodes:
            idx = self.nodes.index(n)
            p = Path(self.nodes[idx:] + self.nodes[:idx], circular=True)
        else:
            p = self.invert()
            idx = p.nodes.index(n)
            p = Path(p.nodes[idx:] + p.nodes[:idx], circular=True)
        return p

    def rotate_to_node(self, node: Node) -> "Path":
        return self.rotate_to(node.id, node.strand)

    def invert(self) -> "Path":
        return Path([n.invert() for n in self.nodes[::-1]], circular=self.circular)

    def __invert__(self) -> "Path":
        return self.invert()

    def __eq__(self, o: object) -> bool:
        return self.nodes == o.nodes

    def __hash__(self) -> int:
        return hash(tuple(self.nodes))

    def __repr__(self) -> str:
        return "_".join([str(n) for n in self.nodes])

    def __len__(self) -> int:
        return len(self.nodes)

    def rename_bids(self, bid_dict: dict) -> "Path":
        return Path(
            [Node(bid_dict[n.id], n.strand) for n in self.nodes], circular=self.circular
        )

    def to_list(self):
        return [n.to_str_id() for n in self.nodes]

    @staticmethod
    def from_list(path_list: list[Node], circular: bool) -> "Path":
        return Path([Node.from_str_id(nid) for nid in path_list])


class Edge:
    """Oriented link between two nodes/paths"""

    def __init__(self, left, right) -> None:
        self.left = left
        self.right = right

    def invert(self) -> "Edge":
        return Edge(self.right.invert(), self.left.invert())

    def __invert__(self) -> "Edge":
        return self.invert()

    def __side_eq__(self, o: object) -> bool:
        return self.left == o.left and self.right == o.right

    def __eq__(self, o: object) -> bool:
        return self.__side_eq__(o) or self.__side_eq__(o.invert())

    def __side_hash__(self) -> int:
        return hash((self.left, self.right))

    def __hash__(self) -> int:
        return self.__side_hash__() ^ self.invert().__side_hash__()

    def __repr__(self) -> str:
        return f"{self.left} <--> {self.right}"

    def __to_str_id(self) -> str:
        return "__".join([self.left.to_str_id(), self.right.to_str_id()])

    def to_str_id(self) -> str:
        A = self.__to_str_id()
        B = self.invert().__to_str_id()
        return A if A < B else B

    @staticmethod
    def from_str_id(t) -> "Edge":
        left, right = t.split("__")
        return Edge(Node.from_str_id(left), Node.from_str_id(right))


def pangraph_to_path_dict(pan):
    """Creates a dictionary isolate -> path objects"""
    res = {}
    for path in pan.paths:
        name = path.name
        B = path.block_ids
        S = path.block_strands
        nodes = [Node(b, s) for b, s in zip(B, S)]
        res[name] = Path(nodes, path.circular)
    return res


def filter_paths(paths, keep_f):
    """Given a filter function, removes nodes that fail the condition from
    the path dictionaries."""
    res = {}
    for iso, path in paths.items():
        filt_path = Path(
            [node for node in path.nodes if keep_f(node.id)],
            path.circular,
        )
        res[iso] = filt_path
    return res


def path_categories(paths):
    """Returns a list of touples, one per non-empty path, with the following info:
    (count, path, [list of isolates])"""
    iso_list = defaultdict(list)
    n_paths = defaultdict(int)
    nodes = {}
    for iso, path in paths.items():
        if len(path.nodes) > 0:
            n_paths[path] += 1
            iso_list[path].append(iso)
            nodes[path] = path.nodes

    # sort by count
    path_cat = [(count, nodes[path], iso_list[path]) for path, count in n_paths.items()]
    path_cat.sort(key=lambda x: x[0], reverse=True)
    return path_cat


def path_edge_count(paths):
    """Count internal edges of paths"""
    ct = Counter()
    for iso, p in paths.items():
        L = len(p.nodes)
        es = []
        for i in range(L - 1):
            e = Edge(p.nodes[i], p.nodes[i + 1])
            es.append(e)
        if p.circular:
            e = Edge(p.nodes[-1], p.nodes[0])
            es.append(e)
        ct.update(es)
    return dict(ct)


def path_block_count(paths):
    """Count internal blocks of paths"""
    ct = Counter()
    for iso, p in paths.items():
        for node in p.nodes:
            ct.update([node.id])
    return dict(ct)


def find_mergers(paths):
    """Create a dictionary source -> sinks of block-ids to be merged"""
    edge_ct = path_edge_count(paths)
    block_ct = path_block_count(paths)

    mergers = {}
    for e, ec in edge_ct.items():
        bl, br = e.left.id, e.right.id
        if (ec == block_ct[bl]) and (ec == block_ct[br]):
            # merge
            if bl in mergers:
                if br in mergers:
                    source = mergers[br]
                    sink = mergers[bl]
                    for k in mergers:
                        if mergers[k] == source:
                            mergers[k] = sink
                else:
                    mergers[br] = mergers[bl]
            elif br in mergers:
                mergers[bl] = mergers[br]
            else:
                mergers[br] = bl
                mergers[bl] = bl

    # add missing blocks that are not in a merger
    for bid in block_ct.keys():
        if bid not in mergers:
            mergers[bid] = bid

    return mergers

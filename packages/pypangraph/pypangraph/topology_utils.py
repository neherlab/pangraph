from collections import Counter, defaultdict


class OrientedBlock:
    """Combination of block id and strandedness"""

    def __init__(self, bid: str, strand: bool) -> None:
        self.id = bid
        self.strand = strand

    def invert(self) -> "OrientedBlock":
        return OrientedBlock(self.id, not self.strand)

    def __invert__(self) -> "OrientedBlock":
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
    def from_str_id(t) -> "OrientedBlock":
        # Block ids are kept as strings internally (they are u64 hashes; see the
        # pypangraph data-model notes), so the id token round-trips as-is.
        bid_str, strand_str = t.split("_")
        return OrientedBlock(bid_str, strand_str == "f")


class Walk:
    """An ordered traversal of oriented blocks (a walk through the block graph)."""

    def __init__(self, nodes=None, circular=None) -> None:
        if nodes is None:
            nodes = []
        self.nodes = nodes
        self.circular = circular

    def add_left(self, oriented_block: OrientedBlock) -> None:
        self.nodes.insert(0, oriented_block)

    def add_right(self, oriented_block: OrientedBlock) -> None:
        self.nodes.append(oriented_block)

    def rotate_to(self, bid: str, strand: bool) -> "Walk":
        if not self.circular:
            raise ValueError("Walk is not circular")
        if bid not in [n.id for n in self.nodes]:
            raise ValueError(f"Block {bid} not in walk")
        n = OrientedBlock(bid, strand)
        if n in self.nodes:
            idx = self.nodes.index(n)
            p = Walk(self.nodes[idx:] + self.nodes[:idx], circular=True)
        else:
            p = self.invert()
            idx = p.nodes.index(n)
            p = Walk(p.nodes[idx:] + p.nodes[:idx], circular=True)
        return p

    def rotate_to_oriented_block(self, oriented_block: OrientedBlock) -> "Walk":
        return self.rotate_to(oriented_block.id, oriented_block.strand)

    def invert(self) -> "Walk":
        return Walk([n.invert() for n in self.nodes[::-1]], circular=self.circular)

    def __invert__(self) -> "Walk":
        return self.invert()

    def __eq__(self, o: object) -> bool:
        return self.nodes == o.nodes

    def __hash__(self) -> int:
        return hash(tuple(self.nodes))

    def __repr__(self) -> str:
        return "_".join([str(n) for n in self.nodes])

    def __len__(self) -> int:
        return len(self.nodes)

    def rename_bids(self, bid_dict: dict) -> "Walk":
        return Walk(
            [OrientedBlock(bid_dict[n.id], n.strand) for n in self.nodes],
            circular=self.circular,
        )


class Edge:
    """Oriented link between two oriented blocks."""

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
        return Edge(OrientedBlock.from_str_id(left), OrientedBlock.from_str_id(right))


def pangraph_to_walks(pan):
    """Creates a dictionary isolate -> Walk objects"""
    res = {}
    for name, path in pan.paths.items():
        B, S = pan.nodes.nodes_to_blocks(path.nodes)
        oriented_blocks = [OrientedBlock(b, s) for b, s in zip(B, S)]
        res[name] = Walk(oriented_blocks, path.circular)
    return res


def filter_walks(walks, keep_f):
    """Given a filter function, removes oriented blocks that fail the condition
    from the walk dictionary."""
    res = {}
    for iso, walk in walks.items():
        filt_walk = Walk(
            [ob for ob in walk.nodes if keep_f(ob.id)],
            walk.circular,
        )
        res[iso] = filt_walk
    return res


def walk_categories(walks):
    """Returns a list of tuples, one per non-empty walk, with the following info:
    (count, walk_nodes, [list of isolates])"""
    iso_list = defaultdict(list)
    n_walks = defaultdict(int)
    nodes = {}
    for iso, walk in walks.items():
        if len(walk.nodes) > 0:
            n_walks[walk] += 1
            iso_list[walk].append(iso)
            nodes[walk] = walk.nodes

    # sort by count
    walk_cat = [(count, nodes[walk], iso_list[walk]) for walk, count in n_walks.items()]
    walk_cat.sort(key=lambda x: x[0], reverse=True)
    return walk_cat


def walk_edge_count(walks):
    """Count internal edges of walks"""
    ct = Counter()
    for iso, w in walks.items():
        L = len(w.nodes)
        es = []
        for i in range(L - 1):
            e = Edge(w.nodes[i], w.nodes[i + 1])
            es.append(e)
        if w.circular:
            e = Edge(w.nodes[-1], w.nodes[0])
            es.append(e)
        ct.update(es)
    return dict(ct)


def walk_block_count(walks):
    """Count occurrences of each block across walks"""
    ct = Counter()
    for iso, w in walks.items():
        for ob in w.nodes:
            ct.update([ob.id])
    return dict(ct)


def find_mergers(walks):
    """Create a dictionary source -> sinks of block-ids to be merged"""
    edge_ct = walk_edge_count(walks)
    block_ct = walk_block_count(walks)

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

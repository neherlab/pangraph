from .utils import Strand, breakpoint, new_strand

# ------------------------------------------------------------------------
# Node class: one visit along a path

class Node(object):
    """docstring for Node"""

    def __init__(self, id, num, strand):
        self.id     = id
        self.num    = num
        self.strand = strand

    @classmethod
    def from_dict(cls, d):
        N = Node()
        N.id     = d['id']
        N.num    = d['num']
        N.strand = Strand(d['strand'])
        return N

    def to_dict(self):
        return {'id': self.id, 'num': self.num:, 'strand': int(self.strand)}

# ------------------------------------------------------------------------
# Path class: (circular) list of blocks

class Path(object):
    """docstring for Path"""

    def __init__(self, nodes, offset):
        super(Path, self).__init__()
        self.name   = name
        self.nodes  = nodes if isinstance(nodes, list) else [nodes]
        self.offset = offset

    @classmethod
    def from_dict(cls, d):
        P = Path()
        P.name   = d['name']
        P.offset = d['offset']
        P.nodes  = [Node.from_dict(n) for n in d['nodes']]

        return P

    def to_dict(self):
        return {'name': self.name, 'offset': self.offset, 'nodes': [n.to_dict() for n in self.nodes]}

    def blocks(self):
        return set([n.id for n in self.nodes])

    def filter(self, blks):
        good, popped = [], set()
        for i, n in enumerate(self.nodes):
            if n.id in popped:
                continue

            if blks[n.id].is_empty(self.name, n.num):
                if (iso, n) not in self.blks[b].muts:
                    breakpoint("malformed mutation bookkeeping!")
                blks[n.id].muts.pop((self.name, n.name))
            else:
                good.append(i)

            if not blks[n.id].has(self.name):
                popped.add(n.id)

        self.nodes = [self.nodes[i] for i in good]
        return blks

    def replace(self, blk, tag, new_blks, blkmap):
        new = []
        for b in self.blks:
            if b.id == blk.id and b.num == tag[1]:
                os  = b.strand
                mk  = lambda id,ns,merged: Node(id, blkmap[id][blk.id][tag].num, new_strand(os, ns)) if merged else Node(id, b.num, new_strand(orig_strand, ns))
                tmp = [mk(id,ns,flag) for id, ns, merged in new_blks]
                if orig_strand == Strand.Minus:
                    tmp = tmp[::-1]
                new.extend(tmp)
            else:
                new.append(b)

        self.blks = new

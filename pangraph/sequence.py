from .utils import Strand, log, breakpoint, new_strand, rev_cmpl

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
        return {'id': self.id, 'num': self.num, 'strand': int(self.strand)}

# ------------------------------------------------------------------------
# Path class: (circular) list of blocks

class Path(object):
    """docstring for Path"""

    def __init__(self, name, nodes, offset):
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

    def sequence(self, blks, verbose=False):
        seq = ""
        for n in self.nodes:
            s = blks[n.id].extract(self.name, n.num, strip_gaps=False, verbose=verbose)
            if n.strand == Strand.Plus:
                seq += s
            else:
                seq += rev_cmpl(s)

        if self.offset != 0:
            seq = seq[self.offset:] + seq[:self.offset]

        return seq

    def rm_empty(self, blks):
        good, popped = [], set()
        for i, n in enumerate(self.nodes):
            if n.id in popped:
                continue

            if blks[n.id].is_empty(self.name, n.num):
                if (self.name, n.num) not in blks[n.id].muts:
                    breakpoint("malformed mutation bookkeeping!")
                blks[n.id].muts.pop((self.name, n.num))
            else:
                good.append(i)

            if not blks[n.id].has(self.name):
                popped.add(n.id)

        self.nodes = [self.nodes[i] for i in good]
        return blks

    def replace(self, blk, tag, new_blks, blk_map, blks):
        new = []
        for b in self.nodes:
            if b.id == blk.id and b.num == tag[1]:
                os  = b.strand
                mk  = lambda id,ns,merged: Node(id, blk_map[id][blk.id][tag][1], new_strand(os, ns)) if merged else Node(id, b.num, new_strand(os, ns))
                tmp = [mk(id,ns,flag) for id, ns, flag in new_blks]
                if os == Strand.Minus:
                    tmp = tmp[::-1]
                new.extend(tmp)
            else:
                new.append(b)

        self.nodes = new

import sys
import numpy as np

from .utils import Strand, log, breakpoint, new_strand, rev_cmpl

# ------------------------------------------------------------------------
# Node class: one visit along a path

class Node(object):
    """docstring for Node"""

    def __init__(self, blk, num, strand):
        self.blk    = blk
        self.num    = num
        self.strand = strand

    def __str__(self):
        return f"({self.blk}, {self.num}, {self.strand})"

    def __repr__(self):
        return str(self)

    @classmethod
    def from_dict(cls, d, blks):
        N = Node()
        N.blk    = blks[d['id']]
        N.num    = d['num']
        N.strand = Strand(d['strand'])
        return N

    def to_dict(self):
        return {'id': self.blk.id, 'num': self.num, 'strand': int(self.strand)}

    def length(self, name):
        return self.blk.len_of(name, self.num)

# ------------------------------------------------------------------------
# Path class: (circular) list of blocks

class Path(object):
    """docstring for Path"""

    def __init__(self, name, nodes, offset):
        super(Path, self).__init__()
        self.name     = name
        self.nodes    = nodes if isinstance(nodes, list) else [nodes]
        self.offset   = offset
        self.position = np.cumsum([0] + [n.length(name) for n in self.nodes])

    def __str__(self):
        return f"{self.name}: {[str(n) for n in self.nodes]}"

    def __repr__(self):
        return str(self)

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
        return set([n.blk for n in self.nodes])

    def sequence(self, verbose=False):
        seq = ""
        for n in self.nodes:
            s = n.blk.extract(self.name, n.num, strip_gaps=False, verbose=verbose)
            if n.strand == Strand.Plus:
                seq += s
            else:
                seq += rev_cmpl(s)

        if self.offset != 0:
            seq = seq[self.offset:] + seq[:self.offset]

        return seq

    def position_of(self, blk):
        for i, n in enumerate(self.nodes):
            if n.blk == blk:
                return i, n.num
        raise ValueError("block not found in path")

    def rm_nil_blks(self):
        good, popped = [], set()
        for i, n in enumerate(self.nodes):
            if n.blk.id in popped:
                continue

            if n.blk.is_empty(self.name, n.num):
                if (self.name, n.num) not in n.blk.muts:
                    breakpoint("malformed mutation bookkeeping!")
                n.blk.muts.pop((self.name, n.num))
            else:
                good.append(i)

            if not n.blk.has(self.name):
                popped.add(n.blk.id)

        self.nodes    = [self.nodes[i] for i in good]
        self.position = np.cumsum([0] + [n.length(self.name) for n in self.nodes])

    def replace(self, blk, tag, new_blks, blk_map):
        new = []
        for n in self.nodes:
            if n.blk.id == blk.id and n.num == tag[1]:
                os  = n.strand
                mk  = lambda b,ns,merged: Node(b, blk_map[b.id][blk.id][tag][1], new_strand(os, ns)) if merged else Node(b, n.num, new_strand(os, ns))
                tmp = [mk(blk,ns,flag) for blk, ns, flag in new_blks]
                if os == Strand.Minus:
                    tmp = tmp[::-1]
                new.extend(tmp)
            else:
                new.append(n)

        self.nodes    = new
        self.position = np.cumsum([0] + [n.length(self.name) for n in self.nodes])

    def __getitem__(self, index):
        if isinstance(index, slice):
            beg = index.start or 0
            end = index.stop or self.position[-1]

            i = np.searchsorted(self.position, beg, side='right')
            j = np.searchsorted(self.position, end, side='right') + 1
            assert i < j, "sorted"
            return [n.blk for n in self.nodes[i:j]]
        elif isinstance(index, int):
            i = np.searchsorted(self.position, index, side='left')
            return self.nodes[i].blk
        else:
            raise ValueError(f"type '{type(index)}' not supported as index")

    def __len__(self):
        return self.position[-1]

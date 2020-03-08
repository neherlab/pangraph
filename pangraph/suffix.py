import io
import pickle
import numpy as np
import suffix_tree

from .utils import Strand, Complement

def base(name):
    return name[:-4]

def direction(name):
    return 1 if name[-3:] == "fwd" else -1

def isfwd(name):
    return name[-3:] == "fwd"

def cat(n, term=True):
    e = ("$", Strand.Null)
    path = n.path
    if path.p < 0:
        return [e]

    if term:
        return [c if isinstance(c, tuple) else e for c in path.S][path.start:path.end]
    else:
        return [c for c in path.S if isinstance(c, tuple)][path.start:path.end]

def plen(n):
    e = ("$", Strand.Null)
    path = n.path
    if path.p < 0:
        return 0
    else:
        return path.end-path.start

def revcmp(seq):
    return [(c[0], Complement(c[1])) for c in seq[::-1]]

def isparent(p, c):
    for ck, cv in c.ident.items():
        contin = False
        for pk, pv in p.ident.items():
            if ck == pk and pv == cv:
                contin = True
                break
        if not contin:
            return False

    return True

def isparent2(p, c):
    P = set(p.ident.items())
    C = set(c.ident.items())

    return C.issubset(P)

def splits(n, *args):
    idxs = [set() for _ in range(len(args))]
    assert not all(idxs)

    for i, arg in enumerate(args):
        for j, ident in enumerate(n.ident):
            if ident[0] == arg:
                idxs[i].add(j)

    if not all(idxs):
        return False

    for i, idx in enumerate(idxs[:-1]):
        for cmpr in idxs[i+1:]:
            if idx == cmpr:
                return False

    return True

def issublist(qry, ref):
    for i in range(len(ref)-len(qry) + 1):
        hit = True
        for j, q in enumerate(qry):
            if q != ref[i+j]:
                hit = False
                break
        if hit:
            return True

    return False

def pprint(db):
    def convert(s):
        if s == Strand.Plus:
            return "+"
        elif s == Strand.Minus:
            return "-"
        elif s == Strand.Null:
            return "/"
        else:
            raise ValueError("Invalid strand sign")

    def flatten(s):
        x, y = zip(*s)
        return "".join(x), "".join(convert(s) for s in y)

    for name, fatstr in db.items():
        S, T = flatten(fatstr)

        print(f"name: {name}")
        print(f"--> str: {S}")
        print(f"--> sgn: {T}")

class Tree(object):
    """
    docstring for suffix tree
    """
    def __init__(self, G):
        if len(G) == 0:
            return

        def mk(lbl, seq, fwd=True):
            if fwd:
                return f"{lbl}_fwd", seq + seq[:-1]
            else:
                seq = revcmp(seq)
                return f"{lbl}_rev", seq + seq[:-1]

        strs = dict([mk(*entry, fwd=True) for entry in G.items()])
        strs.update(dict([mk(*entry, fwd=False) for entry in G.items()]))
        # pprint(strs)

        self.strs = strs
        self.tree = suffix_tree.Tree(strs)
        self.seqs = G

        self.seqlen = { s:len(seq) for s, seq in self.seqs.items() }

        # -----------------------------------------------------------
        # Annotation functions

        def putname(n):
            n.ident = {}
            if n.is_leaf():
                n.ident = {(base(n.str_id), direction(n.str_id)): (n.path.end - n.string_depth()) % self.seqlen[base(n.str_id)] }
            else:
                for c in n.children.values():
                    n.ident.update(c.ident)

        def putlen(n):
            if n.is_leaf():
                n.seql = self.seqlen[base(n.str_id)]
            else:
                n.seql = np.max([c.seql for c in n.children.values()])

        self.tree.root.post_order(putlen)
        self.tree.root.post_order(putname)

    def _getstr_(self, n):
        sfx    = n.string_depth()
        strlen = len(self.strs[n.str_id]) + 1
        s = self.strs[n.str_id][strlen-sfx:strlen]
        if isinstance(s, str):
            s += "$"
        else:
            s.append("$")
        return s

    def __str__(self):
        buf = io.StringIO()
        def tostr(n):
            if n.is_leaf():
                s = self._getstr_(n)
                buf.write(f"{n.ident}\n--> {s}\n")

        self.tree.root.post_order(tostr)
        return buf.getvalue()

    def __repr__(self):
        return str(self)

    @classmethod
    def fromdict(cls, d):
        T = Tree()
        T.strs   = d['strs']
        T.tree   = pickle.loads(d['tree'])
        T.seqs   = d['seqs']
        T.seqlen = d['seqlen']

        return T

    def matches(self, *args):
        if len(args) <= 1:
            return None

        if len(set(args)) == 1:
            return self.seqs[args[0]]

        xlen = min(len(self.seqs[arg]) for arg in args)

        mums = []
        def pushmum(n):
            if plen(n) > xlen:
                return

            if splits(n, *args):
                for i, m in enumerate(mums):
                    if n.suffix_link == m:
                        mums[i] = n
                        return

                    # assert isparent(n,m) == isparent2(n,m), "failure of equivalency"
                    if isparent2(n, m):
                        return

                mums.append(n)

        self.tree.root.post_order(pushmum)

        mums.sort(key = lambda m: plen(m), reverse=True)

        # Remove reversed versions of the same match
        def isrot(mum, other, reverse=False):
            fwd = cat(mum, term=False)
            rev = cat(other, term=False)

            tot = fwd + fwd
            if reverse:
                return issublist(revcmp(rev), tot)
            else:
                return issublist(rev, tot)

        delidx = set([])
        for i, m in enumerate(mums):
            if i in delidx:
                continue

            for j, o in enumerate(mums):
                if j == i or j in delidx:
                    continue
                if isrot(m, o) or isrot(m, o, reverse=True):
                    delidx.add(j)

        for idx in sorted(delidx, reverse=True):
            del mums[idx]

        def unpack(n):
            p = n.path
            return p.S[p.start:p.end]

        return [unpack(m) for m in mums]

    def todict(self):
        return {'strs':   self.strs,
                'tree':   pickle.dumps(self.tree),
                'seqs':   self.seqs,
                'seqlen': self.seqlen}

# ------------------------------------------------------------------------
# Unit tests

def name(n):
    return f"iso{n:02d}"

def putstrand(s):
    if isinstance(s, str):
        return [(a, Strand.Plus) for a in s]
    elif isinstance(s, tuple):
        def convert(b):
            if b == "+":
                return Strand.Plus
            elif b == "-":
                return Strand.Minus
            else:
                raise ValueError(f"unrecognized value: {b}")

        assert len(s) == 2 and len(s[0]) == len(s[1]) and \
                isinstance(s[0], str) and isinstance(s[1], str), "bad input"
        return [(a, convert(b)) for a, b in zip(s[0], s[1])]
    else:
        raise TypeError(f"invalid input: s is type '{type(s)}'")

def label(egs):
    return { name(n):putstrand(eg) for n, eg in enumerate(egs)}

if __name__ == "__main__":
    eg = ["ABCDEF",
         ("FEDCB", "-----"),
         ("CABDEF", "---+++")]
    eg = label(eg)
    T  = Tree(eg)
    M1 = T.matches(name(0), name(1))
    M2 = T.matches(name(0), name(2))
    M3 = T.matches(name(1), name(2))
    print(f"Matches 1:\n{M1}")
    print(f"Matches 2:\n{M2}")
    print(f"Matches 3:\n{M3}")

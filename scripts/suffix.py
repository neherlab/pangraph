import io
import numpy as np
import suffix_tree

PSD  = +1
MSD  = -1
TERM = None

def base(name):
    return name[:-4]

def direction(name):
    return 1 if name[-3:] == "fwd" else -1

def isfwd(name):
    return name[-3:] == "fwd"

def cat(n, term=True):
    e = ("$", 0)
    path = n.path
    if path.p < 0:
        return [e]

    if term:
        return [c if isinstance(c, tuple) else e for c in path.S][path.start:path.end]
    else:
        return [c for c in path.S if isinstance(c, tuple)][path.start:path.end]

def revcmp(seq):
    return [(c[0], -c[1]) for c in seq[::-1]]

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

# sfx = self
class Tree(object):
    """
    docstring for suffix tree
    """
    def __init__(sfx, G):
        def mk(lbl, seq, fwd=True):
            if fwd:
                return f"{lbl}_fwd", seq + seq[:-1]
            else:
                seq = revcmp(seq)
                return f"{lbl}_rev", seq + seq[:-1]

        strs = dict([mk(*entry, fwd=True) for entry in G.items()])
        strs.update(dict([mk(*entry, fwd=False) for entry in G.items()]))

        sfx.strs = strs
        sfx.tree = suffix_tree.Tree(strs)
        sfx.seqs = G

        sfx._annotate_()

    def _annotate_(sfx):
        sfx.seqlen = { s:len(seq) for s, seq in sfx.seqs.items() }
        def putname(n):
            n.ident = {}
            if n.is_leaf():
                n.ident = {(base(n.str_id), direction(n.str_id)): (n.path.end - n.string_depth()) % sfx.seqlen[base(n.str_id)] }
            else:
                for c in n.children.values():
                    n.ident.update(c.ident)

        def putlen(n):
            if n.is_leaf():
                n.seql = sfx.seqlen[base(n.str_id)]
            else:
                n.seql = np.max([c.seql for c in n.children.values()])

        sfx.tree.root.post_order(putlen)
        sfx.tree.root.post_order(putname)

    def _getstr_(sfx, n):
        sfxlen = n.string_depth()
        strlen = len(sfx.strs[n.str_id]) + 1
        s = sfx.strs[n.str_id][strlen-sfxlen:strlen]
        if isinstance(s, str):
            s += "$"
        else:
            s.append("$")
        return s

    def __str__(sfx):
        buf = io.StringIO()
        def tostr(n):
            if n.is_leaf():
                s = sfx._getstr_(n)
                buf.write(f"{n.ident}\n--> {s}\n")

        sfx.tree.root.post_order(tostr)
        return buf.getvalue()

    def __repr__(sfx):
        return str(sfx)

    def matches(sfx, *args):
        if len(args) <= 1:
            return None

        if len(set(args)) == 1:
            return sfx.seqs[args[0]]

        xlen = min(len(sfx.seqs[arg]) for arg in args)

        mums = []
        def pushmum(n):
            if len(cat(n)) > xlen:
                return

            if splits(n, *args):
                for i, m in enumerate(mums):
                    if n.suffix_link == m:
                        mums[i] = n
                        return

                    if isparent(n, m):
                        return

                mums.append(n)

        sfx.tree.root.post_order(pushmum)

        mums.sort(key = lambda m: len(cat(m)), reverse=True)
        # print("\nNon unique maximal matches...")
        # for i, m in enumerate(mums):
        #     print(f"{i}:", m)

        # Remove reversed versions of the same match
        # assert len(mums) % 2 == 0 Not true! Palindromic sequences.
        def isrot(mum, other, reverse = False):
            fwd = cat(mum, term = False)
            rev = cat(other, term = False)

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

        return mums

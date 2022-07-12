import numpy as np

from collections import Counter

from . import pangraph_projector as pgp


def draw_projection(
    proj: pgp.ProjectedGraph,
    ax,
    color_dict,
    exchange=False,
    overlay=True,
    all_dupl=False,
):

    sources = [proj.MPA, proj.MPB]
    if exchange:
        sources = sources[::-1]

    # capture variables
    LsA, LsB = [x.bl_Ls for x in sources]
    pthA, pthB = [x.pth for x in sources]
    commA, commB = [x.comm for x in sources]
    chunkA, chunkB = [x.chunk_id for x in sources]
    did_A, did_B = [x.dupl_id for x in sources]
    sA, sB = [x.s for x in sources]
    colorA, colorB = ["C0", "C1"]
    strA, strB = [proj.strA, proj.strB]
    if exchange:
        colorA, colorB = colorB, colorA
        strA, strB = strB, strA

    # define angular norm
    Ltot = np.sum(LsA)
    free_space = 0.05
    norm = (Ltot * (1.0 + free_space)) / (2 * np.pi)

    # define duplication list
    cA, cB = Counter(pthA), Counter(pthB)
    if all_dupl:
        dupl = set(list(color_dict.keys()))
    else:
        dupl = set([x for x in cA if cA[x] > 1]) | set([x for x in cB if cB[x] > 1])

    # pos dict
    bl_pos = {}

    # draw first path
    x, y = [], []

    t0 = 0
    K = len(LsA)
    de = Ltot * free_space / (norm * K)
    for k in range(K):
        bl = pthA[k]
        l = LsA[k]
        c = color_dict[bl] if bl in dupl else "lightgray"
        r = 1 if commA[k] else 1.03
        dt = l / norm
        xn, yn = draw_arc(ax, r, t0, dt, color=c)
        t0 += dt + de

        x += list(xn)
        y += list(yn)
        # save common blocks
        if commA[k]:
            if not bl in dupl:
                bl_pos[(bl, chunkA[k])] = [xn, yn, sA[k], t0]
            else:
                bl_pos[(bl, chunkA[k], did_A[k])] = [xn, yn, sA[k], t0]

    if overlay:
        x.append(x[0])
        y.append(y[0])
        ax.plot(x, y, color=colorA, alpha=0.1, lw=4, label=strA, zorder=-1)

    # draw second path
    K = len(LsB)
    x, y = [], []
    # starting point: common gene
    k0 = np.argmax(commB == True)
    t0 = 0
    s0 = True
    for k in np.roll(range(K), -k0):

        bl = pthB[k]
        chk = chunkB[k]
        l = LsB[k]
        sb = sB[k]
        if commB[k]:
            if not bl in dupl:
                xn, yn, s, t0 = bl_pos[(bl, chunkB[k])]
            else:
                xn, yn, s, t0 = bl_pos[(bl, chunkB[k], did_B[k])]

            if s != sb:
                xn, yn = xn[::-1], yn[::-1]
                same = False
                # ax.plot(xn, yn, c='red', ls='', marker='.', alpha=0.1)
            else:
                same = True
                # ax.plot(xn, yn, c='red', ls='', marker='.', alpha=0.1)

            x += list(xn)
            y += list(yn)
            r_noise = np.random.rand() * 0.01
        else:
            r = 0.97 + r_noise
            dt = l / norm
            c = color_dict[bl] if bl in dupl else "lightgray"
            if same:
                xn, yn = draw_arc(ax, r, t0, dt, color=c)
                t0 += dt + de

            else:
                xn, yn = draw_arc(ax, r, t0, -dt, color=c)
                t0 -= dt + de

            x += list(xn)
            y += list(yn)

    if overlay:
        x.append(x[0])
        y.append(y[0])
        ax.plot(x, y, color=colorB, alpha=0.1, lw=4, label=strB, zorder=-1)


def draw_arc(ax, r, t0, dt, **kwargs):
    x = r * np.cos([t0, t0 + dt])
    y = r * np.sin([t0, t0 + dt])
    ax.plot(x, y, **kwargs)
    return x, y

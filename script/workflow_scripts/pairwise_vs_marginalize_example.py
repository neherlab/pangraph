import argparse
import numpy as np
import pypangraph as pp
import matplotlib.pyplot as plt


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="""Summary plots for disagreement along a specific pair of strains"""
    )
    parser.add_argument("--pariw", help="input pairwise pangraph", type=str)
    parser.add_argument("--proj", help="input projection pangraph", type=str)
    parser.add_argument("--pdf", help="output pdf plot", type=str)
    return parser.parse_args()


def shared_vector(Ltot, res):
    """Given a set of intervals (beg, end, shared) and the total genome
    length Ltot, builds a boolean vector of length Ltot which is zero
    on the private regions and one on the shared ones."""
    iso = np.zeros(Ltot, dtype=bool)
    for r in res:
        b, e, s = r
        b = (b - 1) % Ltot
        e = e % Ltot
        if not s:
            continue
        if b <= e:
            iso[b:e] = True
        else:
            iso[b:] = True
            iso[:e] = True
    return iso


def build_share_maps(file_pairwise, file_projected):
    """Given a pair of pangraphs builds share-maps for every isolate. These
    are binary vector having the length of the genome and indicating whether
    a region of the genome is shared or not with the other isolate, according
    to either of the graphs."""

    # graphs and pair of isolates
    Gs = {
        "pairwise": pp.Pangraph.load_json(file_pairwise),
        "projected": pp.Pangraph.load_json(file_projected),
    }
    isolates = Gs["pairwise"].strains()

    # map container
    maps = {}

    for iso in isolates:
        for kind in ["pairwise", "projected"]:
            # pick graph and path
            G = Gs[kind]
            p = G.paths[iso]

            # list of keys to access block properties, and block ids
            Keys = zip(p.block_strains, p.block_nums, p.block_strands)
            Bids = p.block_ids

            Ltot, res = 0, []
            for bid, k in zip(Bids, Keys):
                b = G.blocks[bid]
                # shared or not
                shared = len(b.strains()) > 1
                a = b.alignment
                # positions
                beg, end = a.pos[k]
                # add block length
                Ltot += a.block_occurrence_length(k)
                # append shared interval
                res.append((beg, end, shared))

            # build and append map
            maps[(iso, kind)] = shared_vector(Ltot, res)

    return maps


def plot_disagreement(maps, savename):

    isolates = np.unique([k[0] for k in maps])

    fig, axs = plt.subplot_mosaic(
        """
        AAAAAAAAAGG
        BBBBBBBBBGG
        CCCCCCCCCGG
        DDDDDDDDDHH
        EEEEEEEEEHH
        FFFFFFFFFHH
        """,
        figsize=(25, 10),
    )

    for ni, iso in enumerate(isolates):

        # 1) shared vs private for pairwise and projected
        lett = "A" if ni == 0 else "D"
        ax = axs[lett]
        ax.set_title(iso)
        for kind in ["pairwise", "projected"]:
            m = maps[(iso, kind)]
            f = 1 if kind == "pairwise" else -1
            ax.plot(~m * f)
        ax.set_xlim(0, len(m))
        ax.set_yticks([-1, 0, 1])
        ax.set_yticklabels(["private (pairwise)", "shared", "private (projected)"])

        # 2) avree or disagree vs genome position
        lett = "B" if ni == 0 else "E"
        ax = axs[lett]
        diff = np.logical_xor(maps[(iso, "pairwise")], maps[(iso, "projected")])
        ax.plot(diff, color="k")
        ax.set_xlim(0, len(m))
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["agree", "disagree"])

        # 3) cumulative disagreement distribution
        lett = "C" if ni == 0 else "F"
        ax = axs[lett]
        ax.plot(np.cumsum(diff), color="k")
        ax.set_xlim(0, len(m))
        ax.set_xlabel("genome position")
        ax.set_ylabel("cumul. disagree")
        ax.grid(axis="y", alpha=0.3)

        # 4) disagreeement autocorrelation
        lett = "G" if ni == 0 else "H"
        ax = axs[lett]
        corr = []
        deltas = [int(x) for x in np.logspace(0, 4, 50)]
        for l in deltas:
            corr.append(np.mean(diff * np.roll(diff, l)))
        ax.plot(deltas, corr, color="k")
        ax.set_ylim(bottom=0)
        ax.set_xscale("log")
        ax.set_xlabel("$\Delta$")
        ax.set_ylabel("disagreement autocorr.")

    for k in axs:
        axs[k].spines["top"].set_visible(False)
        axs[k].spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # build share/private maps
    maps = build_share_maps(args.pairw, args.proj)

    # summary plots
    plot_disagreement(maps, args.pdf)

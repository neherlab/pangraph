from . import msu_utils as ut
from collections import defaultdict
import numpy as np


def core_paths(pan, L_thr):
    bdf = pan.to_blockstats_df()
    paths = ut.pangraph_to_path_dict(pan)

    def is_core(node_id):
        return (bdf.loc[node_id, "len"] >= L_thr) and bdf.loc[node_id, "core"]

    return ut.filter_paths(paths, is_core)


def flip_msu_to_most_common_orientation(paths):
    orient = defaultdict(int)
    for iso, p in paths.items():
        for n in p.nodes:
            msu_id, strand = n.id, n.strand
            orient[msu_id] += 1 if strand else -1

    # flip all the ones with orient < 0
    for iso, p in paths.items():
        nodes = [n.invert() if orient[n.id] < 0 else n for n in p.nodes]
        paths[iso] = ut.Path(nodes, p.circular)

    return paths


def minimal_synteny_units(pan, L_thr: int, rotate: bool = True):
    c_paths = core_paths(pan, L_thr)
    mergers = ut.find_mergers(c_paths)

    # MSU lengths
    B_len = pan.to_blockstats_df()["len"].to_dict()
    MSU_len = defaultdict(int)
    for bid, msu in mergers.items():
        MSU_len[msu] += B_len[bid]

    # order by length
    MSU_order = sorted(MSU_len, key=MSU_len.get, reverse=True)

    # simplify paths
    MSU_paths = ut.filter_paths(c_paths, lambda x: x in MSU_order)

    # rename MSUs
    MSU_ids = {msu: f"MSU_{i}" for i, msu in enumerate(MSU_order)}
    MSU_len = {MSU_ids[msu]: MSU_len[msu] for msu in MSU_order}
    MSU_paths = {iso: path.rename_bids(MSU_ids) for iso, path in MSU_paths.items()}
    MSU_mergers = {source: MSU_ids[sink] for source, sink in mergers.items()}

    if rotate:
        assert np.all(
            [p.circular for p in MSU_paths.values()]
        ), "Only circular paths can be rotated"
        focal_block = max(MSU_len, key=MSU_len.get)
        focal_strandedness = True
        MSU_paths = {
            iso: p.rotate_to(focal_block, focal_strandedness)
            for iso, p in MSU_paths.items()
        }

    flip_msu_to_most_common_orientation(MSU_paths)

    return MSU_mergers, MSU_paths, MSU_len

import numpy as np


class Path:
    """Python wrapper for pangraph path object. It has attributes:
    - name (str): strain name
    - offset (int): offset of the first block in path?
    - circular (bool): whether the path object is circular

    then there are three attributes that are lists. Each index represent a
    particular block occurrence. These are
    - block_ids (str) : list of block ids
    - block_nums (int) : list of block occurrence (blocks can occurr more than
        once per strain)
    - block_strains (str) : list of strains in which blocks occurr
    - block_strands (bool) : whether the block occurrs on direct or reverse strand.
    - block_positions (int) : list of nucleotide start positions of blocks in the path
    """

    def __init__(self, pan_path):
        self.name = pan_path["name"]
        self.offset = pan_path["offset"]
        self.circular = pan_path["circular"]
        self.block_positions = np.array(pan_path["position"])
        blocks = pan_path["blocks"]
        self.block_ids = np.array([block["id"] for block in blocks])
        self.block_nums = np.array([block["number"] for block in blocks])
        self.block_strains = np.array([block["name"] for block in blocks])
        self.block_strands = np.array([block["strand"] for block in blocks])

    def __len__(self):
        return len(self.block_ids)

    def __str__(self):
        return f"path {self.name}, n. blocks = {len(self.block_ids)}"

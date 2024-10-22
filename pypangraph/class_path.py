from .indexed_collection import IndexedCollection


class Path:
    """Pangraph path object. It has attributes:
    - name (str): strain name
    - circular (bool): whether the path object is circular
    - nodes (list): list of node ids in the path
    - nuc_len (int): total length of the path in base-pairs
    """

    def __init__(self, pan_path):
        self.name = pan_path["name"]
        self.id = pan_path["id"]
        self.circular = pan_path["circular"]
        self.nodes = pan_path["nodes"]
        self.nuc_len = pan_path["tot_len"]

    def __len__(self):
        """Returns the number of nodes in the path"""
        return len(self.block_ids)

    def __str__(self):
        return f"path {self.name}, n. blocks = {len(self.block_ids)}"


class PathCollection(IndexedCollection):
    """Collection of all paths. Inherits from IndexedCollection to allow for
    smart indexing of paths.
    """

    def __init__(self, pan_paths):
        ids = [path["name"] for path in pan_paths.values()]
        items = [Path(path) for path in pan_paths.values()]
        IndexedCollection.__init__(self, ids, items)

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
        return len(self.nodes)

    def __repr__(self):
        return f"path object | name = {self.name}, n. nodes = {len(self.nodes)}, length = {self.nuc_len} bp"

    def __str__(self):
        return f"path object | name = {self.name}, n. nodes = {len(self.nodes)}, length = {self.nuc_len} bp"


class PathCollection(IndexedCollection):
    """Collection of all paths. Inherits from IndexedCollection to allow for
    smart indexing of paths.
    """

    def __init__(self, pan_paths):
        ids = [path["name"] for path in pan_paths.values()]
        # raise an error if there are duplicated path names
        if len(ids) != len(set(ids)):
            raise ValueError("Duplicated path names found in the input json file.")
        items = [Path(path) for path in pan_paths.values()]
        IndexedCollection.__init__(self, ids, items)
        self.idx_to_name = {path.id: path.name for path in items}

import numpy as np
from .class_block import Block
from .class_path import Path
from .class_node import Node


class IndexedCollection:
    """This class is used to implement smart indexing of a list of blocks or paths.
    The class has three elements:
    - ids: ordered list of item ids (string)
    - list: ordered list of items
    - id_to_pos: dictionary mapping ids to their position in the list.

    An element of the class can be indexed in these different ways:
    - through its id (string)
    - through its position on the list
    - though a list of ids
    - through a list of positions
    - though a boolean mask on the list

    This is handled by the __getitem__ function.

    Moreover the object can be cast to iterators, in which case an iterator over
    the list items is returned.

    The object's `len` is the lentgth of the list.
    """

    def __init__(self, ids, items):
        self.ids = np.array(ids)
        self.list = np.array(items)
        self.id_to_pos = {id: n for n, id in enumerate(ids)}

    def __iter__(self):
        return iter(self.list)

    def __len__(self):
        return len(self.list)

    def __getitem__(self, idx):
        # if indexed by block id
        if isinstance(idx, str):
            pos = self.id_to_pos[idx]
            return self.list[pos]

        # if indexed by integer
        if isinstance(idx, (int, np.integer)):
            return self.list[idx]

        # if indexed by list or numpy array
        if isinstance(idx, (list, np.ndarray)):
            # if list is empty return empty list
            if len(idx) == 0:
                return []

            idx0 = idx[0]
            # if the type is integer, return corresponding items
            if isinstance(idx0, (int, np.integer)) and not isinstance(idx0, bool):
                return self.list[idx]

            # if the type is string, return corresponding ids
            if isinstance(idx0, str):
                return [self.__getitem__(idx_i) for idx_i in idx]

            # if the type is bool (a mask)
            if isinstance(idx0, np.bool_):
                return self.list[idx]

        # if no condition is matched, then raise an error
        message = """
        the index object passed does not match any of the allowed types:
        - integer or string
        - list of integers or strings
        - boolean numpy array (mask)
        """
        raise TypeError(message)

    def ids_copy(self):
        return self.ids.copy()


class BlockCollection(IndexedCollection):
    """Collection of all blocks. Inherits from IndexedCollection to allow for
    smart indexing of blocks.
    """

    def __init__(self, pan_blocks):
        ids = [block["id"] for block in pan_blocks]
        items = [Block(block) for block in pan_blocks]
        IndexedCollection.__init__(self, ids, items)


class PathCollection(IndexedCollection):
    """Collection of all paths. Inherits from IndexedCollection to allow for
    smart indexing of paths.
    """

    def __init__(self, pan_paths):
        ids = [path["name"] for path in pan_paths]
        items = [Path(path) for path in pan_paths]
        IndexedCollection.__init__(self, ids, items)

    def to_block_dict(self):
        return {path.name: path.block_ids.copy() for path in self}


class NodeCollection(IndexedCollection):
    """Collection of all nodes. Inherits from IndexedCollection to allow for
    smart indexing of nodes.
    """

    def __init__(self, pan_nodes):
        ids = [node["name"] for node in pan_nodes]
        items = [Node(node) for node in pan_nodes]
        IndexedCollection.__init__(self, ids, items)

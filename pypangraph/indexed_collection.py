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

    The object can be:
    - indexed by id, and the corresponding item is returned.
    - iterated over, yielding the id and the item, as in a dictionary.
    - queried for the list of ids with the `keys` function (like a dictionary).

    The object's `len` is the lentgth of the list.
    """

    def __init__(self, ids, items):
        self.ids = np.array(ids)
        self.list = np.array(items)
        self.id_to_pos = {id: n for n, id in enumerate(ids)}

    def __contains__(self, id_):
        """Returns whether the id is in the collection"""
        return id_ in self.id_to_pos

    def __iter__(self):
        """Returns an iterator over the items, like a dictionary key-value pair"""
        return iter(zip(self.ids, self.list))

    def __len__(self):
        return len(self.list)

    def __getitem__(self, id_):
        """Returns the item corresponding to the id"""
        try:
            pos = self.id_to_pos[id_]
            return self.list[pos]
        except KeyError:
            raise KeyError(
                f"Id {id_} not found in collection {self.__class__.__name__}"
            )

    def keys(self):
        """Returns the list of ids (like a dictionary)"""
        return self.ids.copy()


class BlockCollection(IndexedCollection):
    """Collection of all blocks. Inherits from IndexedCollection to allow for
    smart indexing of blocks.
    """

    def __init__(self, pan_blocks):
        ids = [block["id"] for block in pan_blocks.values()]
        items = [Block(block) for block in pan_blocks.values()]
        IndexedCollection.__init__(self, ids, items)


class PathCollection(IndexedCollection):
    """Collection of all paths. Inherits from IndexedCollection to allow for
    smart indexing of paths.
    """

    def __init__(self, pan_paths):
        ids = [path["name"] for path in pan_paths.values()]
        items = [Path(path) for path in pan_paths.values()]
        IndexedCollection.__init__(self, ids, items)


class NodeCollection(IndexedCollection):
    """Collection of all nodes. Inherits from IndexedCollection to allow for
    smart indexing of nodes.
    """

    def __init__(self, pan_nodes):
        ids = [node["id"] for node in pan_nodes.values()]
        items = [Node(node) for node in pan_nodes.values()]
        IndexedCollection.__init__(self, ids, items)

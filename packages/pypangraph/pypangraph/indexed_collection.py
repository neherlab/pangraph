class IndexedCollection:
    """This class is used to implement smart indexing of a list of blocks or paths.
    The class has three elements:
    - ids: ordered list of item ids (string)
    - list: ordered list of items
    - id_to_pos: dictionary mapping ids to their position in the list.

    The object can be:
    - indexed by id, and the corresponding item is returned.
    - iterated over all of its items like a list.
    - iterated over all the key-value pairs like a dictionary wit the `items` function.
    - queried for the list of ids with the `keys` function (like a dictionary).

    The object's `len` is the lentgth of the list.
    """

    def __init__(self, ids, items):
        self.ids = ids
        self.list = items
        self.id_to_pos = {id_: n for n, id_ in enumerate(ids)}

    def __contains__(self, id_):
        """Returns whether the id is in the collection"""
        return id_ in self.id_to_pos

    def __iter__(self):
        """Returns an iterator over the paths"""
        return iter(self.list)

    def items(self):
        """Returns an iterator over the items, like a dictionary key-value pair"""
        return zip(self.ids, self.list)

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

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

    Ids are stored as strings, but lookups (`__getitem__`, `__contains__`) coerce the
    key to str first, so a collection keyed by numeric ids (e.g. blocks) can be indexed
    with either an int or a str.
    """

    def __init__(self, ids, items):
        self.ids = ids
        self.list = items
        self.id_to_pos = {id_: n for n, id_ in enumerate(ids)}

    def __contains__(self, id_):
        """Returns whether the id is in the collection (accepts int or str)"""
        return str(id_) in self.id_to_pos

    def __iter__(self):
        """Returns an iterator over the paths"""
        return iter(self.list)

    def items(self):
        """Returns an iterator over the items, like a dictionary key-value pair"""
        return zip(self.ids, self.list)

    def __len__(self):
        return len(self.list)

    def __getitem__(self, id_):
        """Returns the item corresponding to the id (accepts int or str)"""
        try:
            pos = self.id_to_pos[str(id_)]
            return self.list[pos]
        except KeyError:
            # `from None`: the original KeyError only carries the key, which is
            # already in this message, so suppress it for a clean traceback.
            raise KeyError(
                f"Id {id_} not found in collection {self.__class__.__name__}"
            ) from None

    def keys(self):
        """Returns the list of ids (like a dictionary)"""
        return self.ids.copy()

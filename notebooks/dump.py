import json
from json import JSONEncoder


class UniversalEncoder(JSONEncoder):
    def default(self, obj):
        if hasattr(obj, '__dict__'):
            return obj.__dict__
        elif hasattr(obj, '__iter__'):
            return list(obj)
        return str(obj)


def dump(obj, fname=None):
    s = json.dumps(obj, cls=UniversalEncoder, indent=2)

    print(s)

    with open(fname or "dump.json", "w") as f:
        f.write(s)

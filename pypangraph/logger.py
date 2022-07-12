# This file implements a simple logger object. This object is to be used in a
# with statement.

from datetime import datetime


class logger:
    """Simple logger class to save the print result to file. Can be used in a
    with statement. The file is created at the beginning of the statement, and
    records are flushed when the statement exits.

    with logger('log.txt', verbose=True) as log:
       log.print('blah')
    """

    def __init__(self, filename, verbose=False):
        self.s = ""
        self.verbose = verbose
        self.filename = filename

    def __enter__(self):
        # initialize clean file
        with open(self.filename, "w") as f:
            now = datetime.now()
            dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
            f.write(f"logger run at {dt_string}\n")
            f.write("---------------------------\n\n")
        return self

    def __exit__(self, type, value, tb):
        self.flush()

    def print(self, *args):
        msg = "  ".join([str(a) for a in args])
        if self.verbose:
            print(msg)
        self.s += str(msg) + "\n"

    def flush(self):
        with open(self.filename, "a") as f:
            f.write(self.s)
        self.s = ""

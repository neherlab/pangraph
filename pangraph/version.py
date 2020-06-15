"""
print the version of pangraph.
"""

from .__version__ import __version__

def register_args(parser):
    pass

def main(args):
    print("pangraph", __version__)
    return 0

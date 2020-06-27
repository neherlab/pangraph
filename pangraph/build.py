"""
build a pangenome alignment from an annotated guide tree
"""
import os, sys
import builtins

from .utils import mkdir
from .tree import Tree

def open(path, *args, **kwargs):
    if path == '-':
        return sys.stdin
    return builtins.open(path, *args, **kwargs)

def register_args(parser):
    parser.add_argument("-d", "--dir",
                        type=str,
                        default=".",
                        help="directory used for output files")
    parser.add_argument("input",
                        type=str,
                        default="-",
                        help="input guide tree [json]")

def main(args):
    '''
    Parameters
    ----------
    args : namespace
        arguments passed in via the command-line from pangraph
    Returns
    -------
    int
        returns 0 for success, 1 for general error
    '''
    with open(args.input, 'r') as input:
        T = Tree.from_json(input)

    tmp = f"{args.dir.rstrip('/')}/tmp"
    mkdir(tmp)
    T.align(tmp)

    return 0

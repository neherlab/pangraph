"""
generates synthetic data from a recombinative Wright-Fisher model
"""

import os

def register_args(parser):
    parser.add_argument("-d", "--dir",
                        type=str,
                        nargs=1,
                        default="data/synth",
                        help="root directory of output files")
    parser.add_argument("-n", "--name",
                        type=str,
                        nargs=1,
                        help="base name for all output files")
    parser.add_argument("-N", "--size",
                        type=int,
                        nargs=1,
                        default=100,
                        help="simulated population size")
    parser.add_argument("-L", "--len",
                        type=int,
                        nargs=1,
                        default=100000,
                        help="average genome length")
    parser.add_argument("-m", "--mu",
                        type=float,
                        nargs=1,
                        default=0,
                        help="average mutation/genome/generation")
    parser.add_argument("-T", "--gen",
                        type=int,
                        nargs=1,
                        default=35,
                        help="number of simulated generations")

def run(args):
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

    return 0

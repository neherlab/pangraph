"""
The top-level pangraph command which dispatches to subcommands.
"""

import argparse
import re
import os
import sys
import importlib

from types import SimpleNamespace

# ------------------------------------------------------------------------
# define subcommands here

subcmd_names = [
    "build",
    "cluster",
    "generate",
    "version",
]

SUB_CMDS = [importlib.import_module('pangraph.' + cmd) for cmd in subcmd_names]

# ------------------------------------------------------------------------
# helper functions

def command_name(cmd):
    """
    Returns a short name for a command module.
    """

    def remove_prefix(prefix, string):
        return re.sub('^' + re.escape(prefix), '', string)

    package     = cmd.__package__
    module_name = cmd.__name__

    return remove_prefix(package, module_name).lstrip(".").replace("_", "-")

def add_default_command(parser):
    """
    Sets the default command to run when none is provided.
    """
    class default_command():
        def main(args):
            parser.print_help()
            return 2

    parser.set_defaults(__command__ = default_command)

def add_version_alias(parser):
    """
    Add --version as a (hidden) alias for the version command.
    It's not uncommon to blindly run a command with --version as the sole
    argument, so its useful to make that Just Work.
    """

    class version_cmd(argparse.Action):
        def __call__(self, *args, **kwargs):
            opts = SimpleNamespace()
            sys.exit( version.main(opts) )

    return parser.add_argument(
        "-v", "--version",
        nargs  = 0,
        help   = argparse.SUPPRESS,
        action = version_cmd
    )

def make_parser():
    parser = argparse.ArgumentParser(
        prog        = "pangraph",
        description = "pangraph: a bioinformatics toolkit for fast multiple genome alignment"
    )

    subparsers = parser.add_subparsers()

    add_default_command(parser)
    add_version_alias(parser)

    for cmd in SUB_CMDS:
        subparser = subparsers.add_parser(
            command_name(cmd),
            help        = (cmd.__doc__).strip().splitlines()[0],
            description = cmd.__doc__)

        subparser.set_defaults(__command__ = cmd)

        cmd.register_args(subparser)

        # use the same formatting class for every command for consistency
        # set here to avoid repeating it in every command's register_parser()
        subparser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    return parser

def main(argv):
    args = make_parser().parse_args(argv)
    return args.__command__.main(args)

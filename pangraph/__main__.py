"""
Stub function and module used as a setuptools entry point.
"""

import pangraph
from sys import argv, exit

# entry point
def main():
    return pangraph.main( argv[1:] )

# if called from python
if __name__ == "__main__":
    exit( main() )

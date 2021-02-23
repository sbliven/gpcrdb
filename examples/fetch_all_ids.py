#!/usr/bin/python
"""
Fetch current list of all GPCRdb structures

@author Spencer Bliven <spencer.bliven@gmail.com>
"""

import sys
import argparse
import logging

# from gpcrdb import


def main(args=None):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "outfile",
        type=argparse.FileType("w"),
        default=sys.stdout,
        nargs="?",
        help="Output file (default: stdout)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Long messages",
        dest="verbose",
        default=False,
        action="store_true",
    )
    args = parser.parse_args(args)

    logging.basicConfig(
        format="%(levelname)s: %(message)s",
        level=logging.DEBUG if args.verbose else logging.WARN,
    )


if __name__ == "__main__":
    main()

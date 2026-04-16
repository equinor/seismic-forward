"""Command-line interface for seismic-forward.
This module provides the command-line interface for running seismic-forward
simulations.
"""

import argparse

import sys
from typing import NoReturn
from importlib.metadata import version, PackageNotFoundError
from .simulation import run_simulation, SeismicForwardError


def get_version() -> str:
    """Get the version of the seismic-forward package.

    Returns:
        str: The version string, or "unknown" if not available.
    """
    try:
        return version("seismic-forward")
    except PackageNotFoundError:
        return "unknown"


def main() -> NoReturn:
    """Main entry point for the CLI.

    Returns:
        NoReturn: The function either exits successfully or with an error code.
    """
    parser = argparse.ArgumentParser(
        prog="seismic_forward",
        description="Seismic Forward Modeling Tool - Generate synthetic seismic from elastic parameters",
    )
    parser.add_argument("modelfile", nargs="?", help="Path to the XML model file")
    parser.add_argument(
        "--version", action="version", version=f"seismic_forward {get_version()}"
    )

    args = parser.parse_args()

    if args.modelfile is None:
        parser.print_help()
        print("\nError: A modelfile must be provided.")
        sys.exit(1)

    try:
        run_simulation(args.modelfile, capture_output=False)
        sys.exit(0)
    except (SeismicForwardError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

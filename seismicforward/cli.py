"""Command-line interface for seismic-forward.
This module provides the command-line interface for running seismic-forward
simulations.
"""
import sys
from typing import NoReturn
from .simulation import run_simulation, SeismicForwardError


def main() -> NoReturn:
    """Main entry point for the CLI.
    
    Returns:
        NoReturn: The function either exits successfully or with an error code.
    """
    if len(sys.argv) != 2:
        print("Error: A modelfile must be provided.")
        print(f"Usage: {sys.argv[0]} modelfile")
        sys.exit(1)

    try:
        run_simulation(sys.argv[1], capture_output=False)
        sys.exit(0)
    except (SeismicForwardError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 

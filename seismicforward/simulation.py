"""Core simulation functionality for seismic-forward.
This module provides the main Python API for running seismic-forward simulations.
It handles the execution of the underlying C++ binary and manages simulation results.
"""
from pathlib import Path
import subprocess
from typing import Union, Dict, Any


class SeismicForwardError(Exception):
    """Base exception for seismic-forward errors."""
    pass


def get_binary_path() -> str:
    """Get the path to the seismic-forward binary.
    
    Returns:
        str: The absolute path to the seismic-forward binary.
        
    Raises:
        FileNotFoundError: If the binary cannot be found in the expected location.
    """
    package_dir = Path(__file__).parent.absolute()
    binary_path = package_dir / "bin" / "seismic_forward"

    if not binary_path.exists():
        raise FileNotFoundError(
            f"Binary not found at {binary_path}. "
            "Please ensure seismic-forward is properly installed."
        )

    return str(binary_path)


def run_simulation(
    model_file: Union[str, Path], 
    capture_output: bool = True
) -> Dict[str, Any]:
    """Run a seismic-forward simulation.
    
    Args:
        model_file: Path to the model file (string or Path object).
        capture_output: Whether to capture and return the simulation output.
            If False, output will be printed directly to stdout.
            
    Returns:
        dict: A dictionary containing simulation results with keys:
            - 'success': bool indicating if simulation completed successfully
            - 'output': str containing stdout (if capture_output=True)
            - 'error': str containing stderr (if capture_output=True)
            
    Raises:
        SeismicForwardError: If the simulation fails or encounters an error.
        FileNotFoundError: If the model file or binary cannot be found.
    """
    model_path = Path(model_file)
    if not model_path.exists():
        raise FileNotFoundError(f"Model file '{model_path}' not found.")

    try:
        binary_path = get_binary_path()
        result = subprocess.run(
            [binary_path, str(model_path)],
            check=False,
            text=True,
            capture_output=capture_output
        )

        output = {
            'success': result.returncode == 0,
            'output': result.stdout if capture_output else None,
            'error': result.stderr if capture_output else None
        }

        if not output['success']:
            raise SeismicForwardError(
                f"Simulation failed with return code {result.returncode}"
            )

        return output

    except Exception as e:
        raise SeismicForwardError(f"Error running simulation: {str(e)}") from e 

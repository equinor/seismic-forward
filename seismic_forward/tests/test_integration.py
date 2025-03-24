import shutil
import gzip
from pathlib import Path
import pytest
from seismic_forward.simulation import run_simulation

def test_run_simulation(tmp_path, monkeypatch):
    # Copy required test files from data directory
    data_dir = Path(__file__).parent.parent.parent / "tests" / "data"
    model_file = data_dir / "modelfile_ps.xml"
    grid_file = data_dir / "Input_grid.grdecl.gz"
    surface_file = data_dir / "Input_top_surface_PS.irap"
    timeshift_file = data_dir / "Input_timeshift_minus50_twt_PS.storm"

    # Copy files to temporary directory
    temp_model_file = tmp_path / "modelfile_ps.xml"
    temp_grid_file = tmp_path / "Input_grid.grdecl"  # Note: removed .gz extension
    temp_surface_file = tmp_path / "Input_top_surface_PS.irap"
    temp_timeshift_file = tmp_path / "Input_timeshift_minus50_twt_PS.storm"
    shutil.copy2(model_file, temp_model_file)
    shutil.copy2(surface_file, temp_surface_file)
    shutil.copy2(timeshift_file, temp_timeshift_file)

    # Copy and decompress the gz file
    with gzip.open(grid_file, 'rb') as f_in:
        with open(temp_grid_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Use monkeypatch to change directory
    monkeypatch.chdir(tmp_path)

    # Test successful simulation
    result = run_simulation(temp_model_file)
    assert result['success'] is True
    assert isinstance(result['output'], str)
    assert isinstance(result['error'], str)

    # Test with non-existent file
    with pytest.raises(FileNotFoundError):
        run_simulation("nonexistent_file.xml")

    # Test with capture_output=False
    result = run_simulation(temp_model_file, capture_output=False)
    assert result['success'] is True
    assert result['output'] is None
    assert result['error'] is None

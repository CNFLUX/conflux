import os
from conflux.config import CONFLUX_DB

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'betaDB', path)

def get_package_data(rel_path):
    """
    Get path to data files packaged with conflux.
    Works with both editable and regular installs.

    Args:
        rel_path: Path relative to package root (e.g., 'data/betaDB/Z_to_element.csv')

    Returns:
        Absolute path to the data file
    """
    return os.path.join(_ROOT, rel_path)

# Load the environment path of the nuclear databases in CONFLUX
__all__ = ['CONFLUX_DB', 'get_package_data']

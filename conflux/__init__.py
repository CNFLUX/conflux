import os
from conflux.config import CONFLUX_DB

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'betaDB', path)

# Load the environment path of the nuclear databases in CONFLUX
__all__ = ['CONFLUX_DB']

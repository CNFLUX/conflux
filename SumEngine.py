# fission product and spectra summation engine

# universal modules
import sys
import numpy as np

# local modules
from BetaEngine import BetaEngine
import FPYEngine import FissionModel

class SumEngine:
    def __init__(self, beta, fission):
        

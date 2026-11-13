import numpy as np
from scipy.interpolate import interp1d

from conflux.bsg.Constants import *

def evaluate_bsg_spectrum(library, parameters, isNu):
    """Evaluate the scalar BSG API on broadcast electron/endpoint energies.

    Array shapes include the energy-bin/sample matrix used for endpoint Monte
    Carlo. Scalar inputs still return a scalar, while arrays retain their shape.
    """
    W, W0 = np.broadcast_arrays(parameters['We'], parameters['W0'])
    args = (parameters['A'], parameters['Z'], parameters['R'], parameters['L'], isNu)
    if W.ndim == 0:
        return library.BSG_beta_spectrum(float(W), float(W0), *args)
    result = np.fromiter(
        (library.BSG_beta_spectrum(float(w), float(w0), *args)
         for w, w0 in zip(W.flat, W0.flat)),
        dtype=float, count=W.size,
    )
    return result.reshape(W.shape)

def stableA(Z):
    """Get estimate for stable mass number as a function of Z"""

    # Parameters from Bethe-Weiszacker formula
    aC = 0.71 # MeV
    aA = 23.7 # MeV
    alpha = aC/4/aA

    return 2*Z*(1+alpha*(2*Z)**(2/3))

def getEltonNuclearRadius(A, natUnits=False):
    """Return nuclear charge radius in fm according to Elton formula

    :param A: mass number
    :param natUnits: return result in units where m_e=hbar=c=1 (Default value = False)

    """
    R = 1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A

    if natUnits:
        """ hbar * c / m_e, in fm"""
        R /= NATURAL_LENGTH/1e-15
    return R

def rebin(x1, x2, y2):
    """Rough rebinning to downsample y2 to grid of x1

    :param x1: array, x values to rebin to
    :param x2: array, initial x values
    :param y2: array, initial y values

    """
    x1Min = np.min(x1)
    x1Max = np.max(x1)
    binSize1 = x1[1] - x1[0]
    binSize2 = x2[1] - x2[0]
    if len(y2.shape) == 1:
        rebinnedY2 = np.zeros(len(x1))
    else:
        rebinnedY2 = np.zeros((len(x1), y2.shape[1]))
    for j in range(len(y2)):
        if (x1Max + binSize1 / 2.0) > x2[j] >= (x1Min - binSize1 / 2.0):
            rebinnedY2[int((x2[j] - (x1Min + binSize1 / 2.0)) / binSize1)] += y2[j] * binSize2 / binSize1
    return rebinnedY2


def downsample(x1, x2, y2, method='linear', fillValue=0., axis=0):
    """Downsample (x2, y2) to (x1, new y)

    :param x1: array, target x range
    :param x2: array, initial x range
    :param y2: array, initial y range
    :param method: interpolation method (Default value = 'linear')
    :param fillValue: default value outside of range (Default value = 0.)

    """
    try:
        f = interp1d(x2, y2, kind=method, bounds_error=False, fill_value=fillValue, axis=axis)
        return np.array(f(x1))
    except ValueError:
        return y2

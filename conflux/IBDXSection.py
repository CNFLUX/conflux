# IBD cross section model

import numpy as np

def ibd_xsection(enu):
    """An empirical IBD cross section function
    Phys. Rev. D 60 053003 (Vogel model)

    :param enu: electron anitneutrino energies array, in MeV
    :return: inverse beta decay cross section in m^2
    :rtype: float"""

    epos = enu - 1.29   # positron energy, MeV
    epos[epos < 0] = 0
    ppos = np.sqrt(epos**2-0.511**2)

    # tau_n = 877.75
    # fr = 1.7152
    # xsection = 2*np.pi/0.511**5/(fr*tau_n)*epos*ppos

    xsection = 0.0952 * epos * ppos * 1e-42
    xsection[np.isnan(xsection)] = 0

    return xsection

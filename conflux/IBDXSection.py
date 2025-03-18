# IBD cross section model

import numpy as np

def ibd_xsection_cm2(Enu):
    """
    IBD interaction cross section in cm^2 / proton.

    From: P. Vogel, J. F. Beacom:
    "An empirical IBD cross section function"
    Phys. Rev. D 60 053003
    Equation 9

    :param Enu: electron anitneutrino energies array, in MeV
    :return: inverse beta decay cross section in cm^2
    :rtype: float
    """

    # positron total energy, MeV $E_e^{(0)} = E_\nu - (m_n - m_p)$
    Ee0 = Enu - 1.29333236
    Ee0[Ee0 <= 0.511] = 0.511  # truncate to pe0 = 0

    # positron momentum
    pe0 = np.sqrt(Ee0**2 - 0.511**2)

    xsection = 0.0952 * Ee0 * pe0 * 1e-42
    xsection[np.isnan(xsection)] = 0

    return xsection

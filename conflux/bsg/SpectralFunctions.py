import numpy as np
from scipy.special import gamma, loggamma, spence
from scipy.special import factorial, factorial2

from conflux.bsg.CoulombFunctions import lambda_k

from conflux.bsg.Constants import *

def phase_space(W, W0, numass = 0, **kwargs):
    """Phase space

    :param W: Electron energy in iunits of me c^2
    :param W0: Electron endpoint energy in units of me c^2

    """
    nue = W0-W-numass
    return np.sqrt(W**2-1)*W*nue*np.sqrt(nue**2-numass**2)

def fermi_function(W, Z, R, **kwargs):
    """Traditional Fermi Function

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param R: Nuclear radius in units of the electron Compton wavelength

    """
    f = 1

    if Z == 0:
        return f

    np.seterr(divide='ignore', invalid='ignore', over='ignore')

    g = np.sqrt(1-(ALPHA*Z)**2)
    p = np.sqrt(W**2-1)
    y = ALPHA*Z*W/p

    #We use the traditional Fermi function, i.e. a prefactor 4 instead of 2(1+gamma)
    #This is consistent with the L0 description below

    f = (4
            *np.power(2*p*R, 2*(g-1))
            *np.exp(np.pi*y)
            *(np.abs(gamma(g+1.j*y))/(gamma(1+2*g)))**2)

    return f

def finite_size_L0(W, Z, R, **kwargs):
    """ Dominant electrostatic finite size correction
    Correction to the traditional Fermi function to use a uniformly charged sphere rather than point charge

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param R: Nuclear radius in units of the electron Compton wavelength

    """
    return (1
            -ALPHA*Z*W*R
            +13/60*(ALPHA*Z)**2
            -ALPHA*Z*R/W)

def finite_size_U_fermi(W, Z, **kwargs):
    """Higher-order electrostatic finite size correction
    Change from uniformly charged sphere to Fermi-type charge distribution

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2

    """
    p = np.sqrt(W**2-1)

    a0 = (-5.6e-5
            -4.94e-5*Z
            +6.23e-8*Z**2)
    a1 = (5.17e-6
            +2.517e-6*Z
            +2e-8*Z**2)
    a2 = (-9.17e-8
            +5.53e-9*Z
            +1.25e-10*Z**2)

    return 1+a0+a1*p+a2*p**2

def sirlin_g(W, W0, **kwargs):
    """Sirlin's g function for order alpha radiative corrections

    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2

    """
    g = 0

    beta = np.sqrt(W**2-1)/W

    np.seterr(divide='ignore', invalid='ignore')
    g = (3*np.log(PROTON_MASS_W)
        -3./4.
        +4./beta*(-1*spence(1-(2*beta/(1+beta))))
        +4*(np.arctanh(beta)/beta-1)*((W0-W)/(3*W)-3/2+np.log(2*(W0-W)))
        +np.arctanh(beta)/beta*(2*(1+beta**2)+(W0-W)**2/(6*W**2)))


    return g

def radiative_correction_o1(W, W0, **kwargs):
    """ Order alpha radiative correction to the beta spectrum shape

    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2

    """
    g = sirlin_g(W, W0)

    return ALPHA/2/np.pi*g

def radiative_correction_o2(W, Z, R, **kwargs):
    """Order alpha^2 Z radiaive correction to the beta spectrum shape

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2

    """
    M = NUCLEON_MASS_W
    L = np.sqrt(10)/R

    gE = 0.5772
    k2 = 0.395

    d1f = (np.log(L/M)
            -k2
            -3/(np.sqrt(10)*np.pi)*L/M
            *(0.5+gE+0.5*np.log(10)+np.log(M/L)))

    d2 = (3/(2*np.sqrt(10)*np.pi)*L/M
            *(1-np.pi/(2*np.sqrt(10))*L/M))

    d3 = (3/(np.sqrt(10)*np.pi)
            *GA
            *GM
            *L/M
            *(gE-1+0.5*np.log(10)+np.log(M/L)+np.pi/(4*np.sqrt(10)*L/M)))

    d01d4 = (np.log(M)
            -2/3*np.log(2*W)
            +35./9
            +np.pi**2/6
            -6*np.log(2))

    return abs(Z)*ALPHA**2*(d1f+d2+d3+d01d4)

def radiative_correction_o3(W, Z, W0, R, **kwargs):
    """Radiative correction of order alpha^3 Z^2 to the beta spectrum shape

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param R: Nuclear charge radius in units of the electron Compton wavelength

    """
    gE = 0.5772
    L = np.sqrt(10)/R

    a = np.pi/3-3/(2*np.pi)
    b = 4/(3*np.pi)*(11/4-gE-np.pi**2/6)
    f = np.log(2*W)-5/6
    g = 0.5*(np.log(R)**2-np.log(2*W)**2)+5/3*np.log(2*R*W)

    return Z**2*ALPHA**3*(a*np.log(L/W)+b*f+4*np.pi/3*g-0.649*np.log(2*W0))

def radiative_correction_L(W0, **kwargs):
    """Resummed order alpha^n radiative corrections

    :param W0: Electron endpoint in units of me c^2

    """
    return 1.026725*(1-2*ALPHA/3/np.pi*np.log(2*W0))**2.25

def radiative_correction(W, Z, W0, R, **kwargs):
    """Total radiative correction up to order alpha^3 Z^2 to the beta spectrum shape

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param R: Nuclear charge radius in units of the electron Compton wavelength

    """
    o1 = radiative_correction_o1(W, W0)
    o2 = radiative_correction_o2(W, Z, R)
    o3 = radiative_correction_o3(W, Z, W0, R)
    L = radiative_correction_L(W0)

    return ((1+o1-ALPHA/2/np.pi*3*np.log(PROTON_MASS_W/W0))
            *(L+o2+o3))

def radiative_correction_neutrino(W, W0, **kwargs):
    """Radiative correction to the (anti)neutrino spectrum to order alpha

    :param Wv: (Anti)neutrino energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2

    """
    r = 1

    p = np.sqrt(W**2-1)
    beta = p/W

    r = (1
            +ALPHA/2/np.pi
            *(3*np.log(PROTON_MASS_W)+23/4
                -8/beta*spence(1-(2*beta)/(1+beta))
                +8*(np.arctanh(beta)/beta-1)*np.log(2*W*beta)
                +4*np.arctanh(beta)/beta*((7+3*beta**2)/8-2*np.arctanh(beta))))
    return r

def recoil_fermi(W, W0, A, **kwargs):
    """Kinematic recoil correction to the beta spectrum shape for a Fermi transition

    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param A: Mass number (protons + neutrons) of the nuclear state

    """
    M = A*NUCLEON_MASS_W
    M2 = M**2

    Vr0 = W0 * W0 / 2. / M2 - 11. / 6. / M2
    Vr1 = W0 / 3. / M2
    Vr2 = 2. / M - 4. * W0 / 3. / M2
    Vr3 = 16. / 3. / M2

    result = (1
            +Vr0
            +Vr1/W
            +Vr2*W
            +Vr3*W**2)

    return result

def recoil_gamow_teller(W, W0, A, **kwargs):
    """Kinematic recoil correction to the beta spectrum shape for a Gamow-Teller transition

    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param A: Mass number (protons + neutrons) of the nuclear state

    """
    M = A*NUCLEON_MASS_W
    M2 = M**2

    Ar0 = -2. * W0 / 3. / M - W0 * W0 / 6. / M2 - 77. / 18. / M2
    Ar1 = -2. / 3. / M + 7. * W0 / 9. / M2
    Ar2 = 10. / 3. / M - 28. * W0 / 9. / M2
    Ar3 = 88. / 9. / M2

    result = (1
            +Ar0
            +Ar1/W
            +Ar2*W
            +Ar3*W**2)

    return result

def recoil_Coulomb_fermi(W, Z, W0, A, **kwargs):
    """Coulomb-recoil correction to the beta spectrum shape for a Fermi transition

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param A: Mass number (protons + neutrons) of the nuclear state

    """
    M = A*NUCLEON_MASS_W
    p = np.sqrt(W**2-1)

    return (1
            -ALPHA*Z*np.pi/M/p
            *(1+(W0-W)/(3*W)))

def recoil_Coulomb_gamow_teller(W, Z, W0, A, **kwargs):
    """Coulomb-recoil correction to the beta spectrum shape for a Gamow-Teller transition

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param A: Mass number (protons + neutrons) of the nuclear state

    """
    beta_type = Z/abs(Z)

    M = A*NUCLEON_MASS_W
    p = np.sqrt(W**2-1)

    return (1
            -beta_type*ALPHA*Z*np.pi/M/p
            *(1-1/3*(W0-W)/(3*W)))

def shape_factor_fermi(W, Z, W0, R, **kwargs):
    """Nuclear shape factor for an allowed Fermi transition

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param R: Nuclear charge radius in units of the electron Compton wavelength

    """
    C0 = (-233/630*(ALPHA*Z)**2
            -1/5*(W0*R)**2
            -6/35*ALPHA*Z*W0*R)
    C1 = (-13/35*ALPHA*Z*R
            +4/15*W0*R**2)
    Cm1 = (2/14*W0*R**2
            +1/70*ALPHA*Z*R)
    C2 = -4/15*R**2

    return 1 + C0 + C1*W + Cm1/W + C2*W**2

def shape_factor_gamow_teller(W, Z, W0, R, A, b, c, d, Lambda, **kwargs):
    """Nuclear shape factor for an allowed Gamow-Teller transition

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param R: Nuclear charge radius in units of the electron Compton wavelength
    :param A: Mass number (protons + neutrons) for the nuclear state
    :param b: Weak magnetism form factor at q2=0
    :param c: Gamow-Teller form factor at q2=0
    :param d: Induced tensor form factor at q2=0
    :param Lambda: Induced pseudoscalar form factor at q2=0

    Form factors are in Holstein notation and that of Hayen et al., RMP 90 (2018) 015008

    """
    M = A*NUCLEON_MASS_W
    bt = Z/abs(Z)

    L = Lambda

    C0 = (-1/5*(W0*R)**2
            +4/9*R**2*(1-L/20)
            +1/3*W0/M/c*(-bt*b+d)
            +2/5*ALPHA*Z/(M*R*c)*(bt*2*b+d)
            +2/35*ALPHA*Z*W0*R*(1-L)
            -233/630*(ALPHA*Z)**2)
    C1 = (bt*4/3*b/M/c
            +4/9*W0*R**2*(1-L/10)
            -4/7*ALPHA*Z*R*(1-L/10))
    Cm1 = (-1/3/M/c*(2*bt*b+d)
            -2/45*W0*R**2*(1-L)
            -ALPHA*Z*R/70)
    C2 = -4/9*R**2*(1-L/10)

    return 1 + C0 + C1*W + Cm1/W + C2*W**2

def shape_factor_unique_forbidden(W, L, W0, Z, R, **kwargs):
    """Unique forbidden shape factor

    :param L: int Spin change
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param Z: Proton number of the final nuclear state

    """
    C = 0.
    pe = np.sqrt(W**2.0 - 1.)
    pnu = W0-W
    L = abs(L)
    for k in range(1, L+1, 1):
        C += lambda_k(W, Z, R, k) * pe**(2*(k-1))*pnu**(2*(L-k))/factorial(2*k-1)/factorial(2*(L-k)+1)
    return C

def atomic_screening(W, Z, R, l, **kwargs):
    """Screening correction due to atomic electrons in the final state

    :param Z: Proton number of the final nuclear state
    :param W: Elecron energy in units of me c^2
    :param R: Nuclear charge radius in units of the electron Compton wavelength
    :param l: Shift in electric potential at the origin due to atomic electrons

    """
    beta_type = Z/abs(Z)

    S = 1
    X = 1

    p = np.sqrt(W**2-1)
    Wt = W - beta_type * 0.5 * ALPHA * (abs(Z)-beta_type)*l

    pt = (0.5*p
            +0.5*np.sqrt(p**2-beta_type*2*ALPHA*abs(Z)*Wt*l+0.j))

    y = ALPHA*Z*W/p
    yt = ALPHA*Z*Wt/pt
    g = np.sqrt(1-(ALPHA*Z)**2)

    S = (np.abs(gamma(g+1.j*yt)/gamma(g+1.j*y))**2
            *np.abs(np.exp(loggamma(g+2.j*pt/l)-loggamma(1+2.j*p/l)))**2
            #*np.abs(gamma(g+1.j*2*pt/l))**2/np.abs(gamma(1+1.j*2*p/l))**2
            *np.exp(-np.pi*y)
            *np.power(2*p/l, 2*(1-g)))

    X = (1/(1+0.25*(l/p)**2)
            *(1+1/8*(Wt+g)/Wt*(l/p)**2
                +0.5*g**2/(1+np.sqrt(1-ALPHA*Z*l/(W+1)))**2
                *(W-1)/Wt*(l/p)**2
                *(1-1/8*(1-g)/g*(l/p)**2)))

    return X*S


def atomic_mismatch(W, Z, W0, A, **kwargs):
    """Correction due to non-orthogonality of initial and final electronic states, resulting in shake-up and shake-off

    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    :param A: Nuclear mass number (protons + neutrons)

    """
    beta_type = Z/abs(Z)
    dBdZ2 = (44.200 * np.power(Z - beta_type, 0.41)
            +2.3196e-7 * np.power(Z - beta_type, 4.45))
    K = -0.872 + 1.270 * np.power(abs(Z), 0.097) + 9.062e-11 * np.power(abs(Z), 4.5)
    beta = np.sqrt(W**2-1)/W
    l = 1.83E-3 * K * Z / beta
    M = A * NUCLEON_MASS_W
    vR = np.sqrt(1 - M * M / (M * M + (W0 * W0 - 1) / 4.))
    psi2 = 1 + 2 * ALPHA / beta * (np.arctan(1 / l) - l / 2 / (1 + l * l))
    C0 = -ALPHA**3*Z / beta * l / (1 + l * l) / psi2 #CHECK ALPHA power
    C1 = (2*ALPHA**2*Z * vR / beta
            *((0.5 + l * l) / (1 + l * l) - l * np.arctan(1 / l)) / psi2)

    return 1 - 2 / (W0 - W) * (0.5 * dBdZ2 + 2 * (C0 + C1))

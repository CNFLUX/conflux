# Copyright 2025 Lawrence Livermore National Security, LLC. See the top-level NOTICE file for details.
# Author: Xianyi Zhang

# SPDX-License-Identifier: MIT

"""Public modules."""
import numpy as np
import xml.etree.ElementTree as ET
from copy import deepcopy
from tqdm import tqdm
from collections import namedtuple

"""CONFLUX modules."""
from conflux.config import CONFLUX_DB
from conflux.Basic import Spectrum, integrate_trapezoid
from conflux.bsg.Constants import ELECTRON_MASS_MEV, NATURAL_LENGTH
from conflux.bsg.SpectralFunctions import (phase_space, 
                                           fermi_function, 
                                           finite_size_L0, 
                                           finite_size_L0_simple,
                                           recoil_gamow_teller, 
                                           radiative_correction, 
                                           radiative_correction_neutrino,
                                           recoil_Coulomb_gamow_teller,                               
                                           atomic_screening,
                                           shape_factor_gamow_teller,
                                           shape_factor_unique_forbidden)
from conflux.bsg.FiniteSize import getL0Constants
from conflux.bsg.Functions import getEltonNuclearRadius
from conflux.bsg.Screening import screening_potential

import matplotlib.pyplot as plt
from ctypes import *

libBSG = None

def init_libBSG(libname = "libBSG.so"):
    global libBSG
    libBSG = cdll.LoadLibrary(libname)
    libBSG.BSG_beta_spectrum.argtypes = [c_double, c_double, c_int, c_int, c_double, c_int, c_bool]
    libBSG.BSG_beta_spectrum.restype = c_double
    print("Calculations using BSG library", libname)

#########################################
# Final neutrino and antineutrino spectra

# Electron or neutrino spectrum
def _e_nu_spectrum(p, isNu, numass = 0):
    """
    Calculate the electron or beta spectrum from theory as a function of energy.

    :param W: electron total energy in electron mass units, (m_e + KE)/m_e
    :type W: float
    :param p: A dictionary containing parameters to be used in the
        calculation of the beta spectrum from theory.
        Parameters: {
                    'We': electron total energy in electron mass units, (m_e + KE)/m_e
                    'Z': Z,
                    'A': A,
                    'W0': Spectrum endpoint W
                    'R': nuclear radius,
                    'L': forbiddenness transition,
                    'c': 1.0,
                    'b': bAc*A,
                    'd': 0.0,
                    'Lambda': 0.0,
                    'l': screening_potential(Z)
                }
    :type p: dictionary
    :param isNu: spectrum corrections for neutrino or electron?
    :type isNu: bool
    :param numass: Neutrino mass parameter, defaults to 0
    :type numass: float, optional
    :return: beta spectrum amplitude at given beta energy
    :rtype: float
    """

    np.seterr(divide='ignore', invalid='ignore', over='ignore')

    if libBSG is not None:
        # determine which inputs are arrays
        array_n = 0
        isarray = {k: hasattr(p[k], '__iter__') for k in {'We', 'W0'}}
        for k,b in isarray.items():
            if b:
                n = len(p[k])
                assert array_n == 0 or n == array_n
                array_n = n

        if not array_n:
            res = libBSG.BSG_beta_spectrum(p["We"], p["W0"], p["A"], p["Z"], p["R"], p['L'], isNu)
        else:
            res = np.zeros(array_n)
            for i in range(array_n):
                res[i] = libBSG.BSG_beta_spectrum(
                    p['We'][i] if isarray['We'] else p['We'],
                    p['W0'][i] if isarray['W0'] else p['W0'],
                    p['A'], p['Z'], p['R'], p['L'], isNu)

        return res

    result = (
        phase_space(p["We"], numass=numass, **p)
        * fermi_function(p["We"], **p)
        * finite_size_L0(p["We"], L0Const = getL0Constants(p['Z']), **p)
    )

    if p['Z']:
        result *= (
            recoil_gamow_teller(p["We"], **p)
            * recoil_Coulomb_gamow_teller(p["We"], **p)
            * atomic_screening(p["We"], **p)
        )

        if p['L'] == 0:
            result *= shape_factor_gamow_teller(p["We"], **p)
        else:
            result *= shape_factor_unique_forbidden(p["We"], **p)

    if isNu: result *= radiative_correction_neutrino(p["We"], **p)
    else:    result *= radiative_correction(p["We"], **p)

    return np.nan_to_num(result, nan=0.0)


def electron(ebeta, p, numass=0):
    """
    Calculate the beta spectrum from theory as a function of energy.

    :param ebeta: The energy of the beta. Unit: MeV
    :type ebeta: float
    :param p: A dictionary containing parameters to be used in the
        calculation of the beta spectrum from theory.
        Parameters: {
                    'Z': Z,
                    'A': A,
                    'W0': Spectrum endpoint electron total energy (m_e + KE)/m_e
                    'R': nuclear radius,
                    'L': forbiddenness transition,
                    'c': 1.0,
                    'b': bAc*A,
                    'd': 0.0,
                    'Lambda': 0.0,
                    'l': screening_potential(Z)
                }
    :type p: dictionary
    :param numass: Neutrino mass parameter, defaults to 0
    :type numass: float, optional
    :return: beta spectrum amplitude at given beta energy
    :rtype: float
    """

    p['We'] = ebeta/ELECTRON_MASS_MEV + 1
    return _e_nu_spectrum(p, False, numass)


def neutrino(enu, p, numass=0):
    """
    Calculate the neutrino spectrum from theory as a function of energy.
    
    :param enu: The energy of the electron antineutrino. Unit: MeV
    :type ebeta: float 
    :param p: A dictionary containing parameters to be used in the 
        calculation of the beta spectrum from theory. 
        Parameters: { 
                    'Z': Z,
                    'A': A,
                    'W0': Spectrum endpoint electron total energy (m_e + KE)/m_e
                    'R': nuclear radius,
                    'L': forbiddenness transition,
                    'c': 1.0,
                    'b': bAc*A,
                    'd': 0.0,
                    'Lambda': 0.0,
                    'l': screening_potential(Z)
                }
    :type p: dictionary
    :param numass: Neutrino mass parameter, defaults to 0
    :type numass: float, optional
    :return: beta spectrum amplitude at given beta energy 
    :rtype: float
    
    """

    p['We'] = np.clip(p['W0'] - enu/ELECTRON_MASS_MEV, 1., None)
    return _e_nu_spectrum(p, True, numass)



#####################################
#####################################
#####################################


class BetaBranch:
    """One beta decay spectrum branch, with precalculated spectrum shape for interpolation"""

    Z: int
    """Atomic number of the mother isotope"""
    A: int
    """The Atomic mass of the mother isotope"""
    I: int
    """Isomeric excited state enumerator of the mother isotope"""

    E0: float
    """End-point energy of this branch, unique for each individual beta-unstable isotope"""
    sigma_E0: float
    """Uncertainty of the end-point energy"""

    """The uncertainty of the branching fraction."""
    forbiddenness: int

    def __init__(self, Z: int, A: int, I: int, E0: float, sigma_E0: float, forbiddenness: int, bAc: float = 4.7):

        self.Z = Z
        self.A = A
        self.I = I

        self.E0 = E0
        self.sigma_E0 = sigma_E0

        self.forbiddenness: int = forbiddenness

        self.bAc = bAc

    def GetParameters(self):
        """Get calculation parameters dict"""
        return {
            'Z': self.Z,
            'A': self.A,
            'W0': self.E0/ELECTRON_MASS_MEV + 1,
            'R': getEltonNuclearRadius(self.A) * 1e-15 / NATURAL_LENGTH,
            'L': self.forbiddenness,
            'c': 1.0,
            'b': self.bAc * self.A,
            'd': 0.0,
            'Lambda': 0.0,
            'l': screening_potential(self.Z)
        }


    def BetaSpectrum(self, e, nu_spectrum=False):
        """
        Calculate the unnormalized beta/neutrino spectral shape as a function of energy.

        :param e: the kinetic energy of beta/neutrino (MeV)
        :type e: float
        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to False
        :type nu_spectrum: bool, optional
        :param numass: neutrino mass, defaults to 0
        :type numass: float, optional
        :return: the spectrum amplitude at specific energy
        :rtype: float

        """

        params = self.GetParameters()

        if self.custom_func is not None:
            function = lambda e: ((1-self.mixing)*self.custom_func(e, params)
                                  +self.mixing*self.custom_func(e, params, numass=self.numass))
        elif nu_spectrum:
            function = lambda e: ((1-self.mixing)*neutrino(e, params)
                                  +self.mixing*neutrino(e, params, numass=self.numass))
        else:
            function = lambda e: ((1-self.mixing)*electron(e, params)
                                  +self.mixing*electron(e, params, numass=self.numass))

        return function(e)

    def CalcNormalizedSpectrum(self, nu_spectrum=False, xbins=None, uncertainties: bool = True):
        """
        Calculated unit-normalized beta/neutrino spectrum shape

        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to False
        :type nu_spectrum: bool, optional
        """

        self.xbins = xbins
        if not self.xbins:
            npts = 50
            self.xbins = (np.arange(0, npts)/(npts-1))*self.E0

        self.spectrum = self.BetaSpectrum(self.xbins, nu_spectrum)
        if uncertainties: self.uncertainty = np.array([self.SpectEndpointUncertMC(x, nu_spectrum) for x in self.xbins])

        # Normalization spectrum (calculate full range if not already covered)
        if self.xbins[0] <= 0 and self.xbins[-1] >= self.E0:
            norm = integrate_trapezoid(self.xbins, self.spectrum)
        else:
            nrmx = np.arange(0., self.E0, self.E0/50.)
            nrmspec = self.BetaSpectrum(nrmx, nu_spectrum)
            norm = integrate_trapezoid(nrmx, nrmspec)

        if not norm > 0:
            print("Warning: failed spectrum normalization", norm)
            self.Display()
            self.spectrum = np.zeros(len(self.xbins))
        else:
            self.spectrum /= norm
            if uncertainties: self.uncertainty /= norm

    def __call__(self, E):
        """Interpolate normalized beta spectrum at specified point"""
        return np.interp(E, self.xbins, self.spectrum)

    # display the vital info of branch
    def Display(self):
        """Display vital branch info."""
        print(f"  Branch Z{self.Z}A{self.A}I{self.I} E0 = {self.E0} +\- {self.sigma_E0} MeV, L = {self.forbiddenness}")
        
    # set correlation between this branch and the other branch
    def SetCovariance(self, otherBranch, correlation):
        """
        Calculate the correlation and covariances between this branch and the other branch.
        
        :param otherBranch: The other beta branch whose correlation and covariance you want to calculate
        :type otherBranch: :class:`conflux.BetaEngine.BetaBranch`
        :param correlation: The correlation between the two branches
        :type correlation: float
        
        """
        #If the endpoint energy between the two branches is the same, do nothing
        if self.E0 == otherBranch.E0: return

        #Otherwise, set the correlation of the other branch at that endpoint energy to the correlation of this branch
        self.corr[otherBranch.E0] = correlation
        #Also, calculate the covariance of the other branch from the fractoinal uncertainties of the two branches and their correlation
        self.cov[otherBranch.E0] = (self.sigma_frac * correlation
                                    * otherBranch.sigma_frac)

    def SpectEndpointUncertMC(self, e, nu_spectrum=False, samples = 30):
        """
        Calculate the beta/neutrino spectral shape uncertainty due to enpoint energy uncertainty sigma_E0
        as a function of energy using Monte Carlo sampling.
        
        :param e: The energy of beta/neutrino (MeV).
        :type e: float
        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to False
        :type nu_spectrum: bool, optional
        :param samples: Set the size of MC samples, defaults to 30
        :type samples: int, optional
        :return: The uncertainty of spectrum amplitude at the input energy
        :rtype: float
        """

        if self.sigma_E0 == 0: return 0.

        E0range = np.random.normal(self.E0, self.sigma_E0, samples)
        newParameters = self.GetParameters()
        newParameters['W0'] = E0range/ELECTRON_MASS_MEV + 1

        if (nu_spectrum == True):
            function = lambda e: ((1-self.mixing)*neutrino(e, newParameters)
                                  +self.mixing*neutrino(e, newParameters, numass=self.numass))
        else:
            function = lambda e: ((1-self.mixing)*electron(e, newParameters)
                                  +self.mixing*electron(e, newParameters, numass=self.numass))
            
        if (self.custom_func!=None):
            function = lambda e: ((1-self.mixing)*self.custom_func(e, newParameters)
                                  +self.mixing*self.custom_func(e, newParameters, numass=self.numass))

        return np.std(function(e))




class DecayBranching:
    """Beta decay branches summary for one parent isotope"""

    Z: int
    """The Atomic number of the mother isotope"""
    A: int
    """The Atomic mass of the mother isotope"""
    I: int
    """The isomeric state this decay mode"""
    HL: float
    """Half life of this isotope (s)"""
    Q: float
    """Q value of the decay (MeV) to daughter ground state"""
    name: str
    """Name of the isotope"""

    missing: bool = False
    """Uncertainty of the end-point energy"""
    branches: list
    """[(frac, dfrac, E0), ...] branching fraction for each beta decay energy E0, sorted ascending by frac"""

    def __init__(self, Z, A, I, Q, HL, name):
        """Constructor method."""
        
        self.HL = HL # top-of-chain isotope half-life (seconds)

        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q

        self.name = name
        self.missing = False
        self.branches = []

        #self.corr = {E0: 1}  # correlation with other branches of the same isotope
        #self.cov = {E0: self.sigma_frac**2} # Set the covariance diagonal element to the square of the branch fraction uncertainty

    def AddBranch(self, fraction, sigma_frac, betaBranch):
        self.branches.append((fraction, sigma_frac, betaBranch))

    def SynthesizeMissingBranches(self, missingBranch):
        """Create sythetic branches"""
        self.missing = True
        self.branches = []
        if self.Q <= 0: missingBranch = 1
        for i in range(missingBranch):
            misE0 = (i+1)*max(0., self.Q)/missingBranch
            self.AddBranch(1./missingBranch, 0., BetaBranch(self.Z, self.A, self.I, misE0, 0., 0))

    def Display(self):
        """Display vital isotope property and branch information."""
        print('Beta isotope: '+self.name+', ZAI = %i %i %i'%(self.Z, self.A, self.I) + ', Q = '
            +str(self.Q)+" MeV, " +str(len(self.branches))+" branches")
        for frac, dfrac, branch in self.branches:
            print("*  %g +- %g"%(frac, dfrac))
            branch.Display()

    def decay_fraction(self, begin, end):
        """
        Calculate the proportion of isotope decayed in the given time window.

        :param begin: the begining of the window (s)
        :type begin: float
        :param end: the end of the window (s)
        :type end: float
        :return: The fraction of isotope decayed (from 0 to 1)
        :rtype: float
        """
        return 2**(-begin/self.HL) - 2**(-end/self.HL)

    def CalcCovariance(self, GSF=True):
        """
        Calculate the covaraince matrix for all the beta branches of this isotope.
        
        :param GSF: Determine whether to consider the anti-correlation between the Ground-State-decay branching Fraction (GSF) and other branches, defaults to True
        :type GSF: bool, optional

        """
        assert False # TODO

        MaxBranch = self.MaxBranch()
        # obtain the ground state branch info
        totalFrac = 0
        GSFrac = MaxBranch.frac
        for E0, branch in self.branches.items():
            totalFrac += branch.frac
        restFrac = totalFrac-GSFrac

        if GSF == False or MaxBranch.E0 != self.Q or restFrac == 0:
            for i, branchi in self.branches.items():
                for j, branchj in self.branches.items():
                    branchi.SetCovariance(branchj, 0)
            return

        # calculate the covariance matrix with ground state anticorrelation
        for i, branchi in self.branches.items():
            for j, branchj in self.branches.items():
                if i == self.Q:
                    correlation = -1*branchj.frac/restFrac
                    branchi.SetCovariance(branchj, correlation)
                    branchj.SetCovariance(branchi, correlation)
                elif j != self.Q:
                    branchi.SetCovariance(branchj, 0)

    def SumSpectra(self, Es, nu_spectrum=True, branchErange=(0., 1000)):
        """
        Sum the beta/antineutrino spectra of all branches of this isotope at energies Es
        
        :param Es: energy points (MeV) to evaluate
        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to True
        :type nu_spectrum: bool, optional
        :param branchErange: set the lower and upper limit of the endpoint energy of branches to be summed, defaults to (0, 1000.)
        :type branchErange: two-element list, optional
        """

        self.nbin = len(Es)
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)
        self.spectUnc = np.zeros(self.nbin) # theoretical uncertainty
        self.branchUnc = np.zeros(self.nbin)
        totalUnc = np.zeros(self.nbin)

        # calculate the total uncertaintty and append spectra
        for E0i, branchi in self.branches.items():
            if E0i < branchErange[0] or E0i > branchErange[1]:
                continue

            branchi.CalcNormalizedSpectrum(nu_spectrum)
            si = branchi(Es)
            fi = branchi.frac

            self.spectrum += si * fi

            continue # TODO uncertainty

            self.spectUnc += branchi.uncertainty*fi

            for E0j, branchj in self.branches.items():
                sj = branchj.spectrum
                # fj = branchi.frac
                # dj = branchj.sigma_frac
                cov_bij = branchi.cov[E0j]
                sigma_bij = si * cov_bij * sj
                self.branchUnc += sigma_bij
                if E0i == E0j:
                    totalUnc += (branchi.uncertainty * fi)**2 + sigma_bij
                else:
                    totalUnc += sigma_bij

        self.branchUnc = np.sqrt(self.branchUnc)
        self.uncertainty = np.sqrt(totalUnc)

    def decay_fraction(self, t, HLs):
        """
        Calculate the fraction of decays completed at a generation in a decay chain.
        The list HLs contains all half lives of the this isotope and all its antecedents.
        
        :param t: time for which to calculate decay fraction (in same units as HLs)
        :type t: float
        :param HLs: half-lives of all decays in chain (in same units as t)
        :type HLs: list(int)
        :return: fraction of decays completed at time t (0 for t=0; 1 as t -> infinity)
        :rtype: float
        """

        rate = 1
        for HL in HLs: rate *= 1 - 2**(-t/HL)
        return rate
    


# BetaEngine tallys beta branches in the betaDB and calculate theoretical beta 
# spectra of all tallied branches
# if inputlist is not given, load the entire betaDB from the default betaDB
class BetaEngine:
    """
    BetaEngine tallys beta branches in the betaDB and calculate theoretical beta spectra of all tallied branches.
    """

    inputlist: list[int] = None
    """A list of isotopes. If the inputlist is not given, load the entire betaDB from the default betaDB."""
    isots: dict
    """A dictionary of isotopes. isots contain keys as the ZAI number of the isotope and values being :class:`conflux.BetaEngine.BetaIstp`"""
    targetDB: str = CONFLUX_DB + "/betaDB/ENSDFbetaDB2.xml"

    """The file name of beta decay data base, defaults to CONFLUX_DB+`/betaDB/ENSDFbetaDB2.xml'"""
    custom_func: callable = None
    """Customized beta function, needs to be structued similarly as :meth:`conflux.BetaEngine.neutrino` or :meth:`conflux.BetaEngine.electron`, defaults to None"""
    numass: float = 0
    """To calculate the spectrum with non-zero neutrino mass, give the neutrino mass in the MeV unit."""
    mixing: float = 0
    """To calcualte spectrum with non-zero neutrino mass, provide the mixing (0-1) of the neutrino mass state."""
    
    
    def __init__(self, inputlist=None,  targetDB=CONFLUX_DB + "/betaDB/ENSDFbetaDB2.xml"):
        """Constructor, loading from datase values"""
        self.inputlist = inputlist
        self.isots = { } # {ZAI: DecayBranching}
        self.LoadBetaDB(targetDB)
        
    def LoadBetaDB(self, targetDB: str = CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml", missingBranch: int = 3) -> None:
        """
        Load default or input betaDB to obtain beta decay informtion. A customed DB must follow the same format as the default DB.
        
        :param targetDB: The file name of beta decay data base, defaults to CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml"
        :param missingBranch: Determine how many branches there are in a missing isotope, defaults to 3
        """

        # empty the existing isotope list before filling loading the DB
        self.isots = {}
        
        if self.inputlist is not None:
            print("Loading all beta data from the default betaDB...")
        print("Searching DB: "+targetDB+"...")
        print("Loading spectra of beta branches...")

        for isotope in ET.parse(targetDB).getroot():

            ZAI = int(isotope.attrib['isotope'])

            # skip isotopes not in inputlist (if provided)
            if self.inputlist is not None and ZAI not in self.inputList: continue

            Z = int(ZAI/10000)
            A = int(ZAI%10000/10)
            I = int(ZAI%10)

            Q = float(isotope.attrib['Q'])
            HL = float(isotope.attrib['HL'])
            name = isotope.attrib['name']

            betaIstp = DecayBranching(Z, A, I, Q, HL, name)
            if not Q > 0:
                print("Warning: negative-Q branch")
                betaIstp.Display()

            # some isotopes contain summed branch fraction greater than 1;
            # prepare to normalize
            fracsum = sum([float(branch.attrib['fraction']) for branch in isotope])

            # actually assign values from database to branches
            for branch in isotope:
                E0 = float(branch.attrib['end_point_E'])
                if not E0 > 0: continue

                sigma_E0 = float(branch.attrib['sigma_E0'])
                spin_par_changes = [n.strip() for n in branch.attrib['dJpi'].split(',')]

                # converting spin change to forbidden rransition types
                # if spin_par_change contain '-', it means parity changed
                ftypes = [['0', '1'], ['0-', '1-', '2-'], ['2', '3'],
                            ['3-', '4-'], ['4', '5']]
                firstftypes = [-10, -11, 10]
                forbiddenness = 6
                for i in range(len(ftypes)):
                    for j in range(len(spin_par_changes)):
                        if (spin_par_changes[j] in ftypes[i] and
                            i < abs(forbiddenness)):
                            forbiddenness = -i
                            if i == 1:
                                # TODO need to add specific type for the first order non-unique transition
                                forbiddenness = firstftypes[i]
                            elif spin_par_changes[j] == ftypes[i][-1]:
                                forbiddenness = i

                # assign fraction values to branches
                # normalize if greater than one
                fraction = float(branch.attrib['fraction'])
                sigma_frac = float(branch.attrib['sigma_frac'])
                if fracsum > 1: fraction /= fracsum

                betaBranch = BetaBranch(Z, A, I, E0, sigma_E0, forbiddenness)
                betaIstp.AddBranch(fraction, sigma_frac, betaBranch)

            # if same ZAI reappears, prefer one with more branches...
            if ZAI in self.isots:
                # print(f"\nWarning: Conflicting Duplicate ZAI {ZAI}")
                # self.isots[ZAI].Display()
                # betaIstp.Display()
                if len(self.isots[ZAI].branches) > len(betaIstp.branches): continue

            self.isots[ZAI] = betaIstp

        for i in self.isots.values():
            i.branches.sort(key=(lambda x: x[0]))
            if not len(i.branches):
                #print("Synthesizing branches for"); i.Display()
                i.SynthesizeMissingBranches(missingBranch)

    def CalcBetaSpectra(self, nu_spectrum: bool = False, uncertainties: bool = True, silent: bool = False):
        """Precalculate spectra shapes for decay branches"""

        for i in tqdm(self.isots.values(),
                      desc="Calculating beta/neutrino spectra of "+str(len(self.isots))+ " isotopes",
                      disable=silent):
            for f,df,b in i.branches:
                b.numass = self.numass
                b.custom_func = self.custom_func
                b.mixing = self.mixing
                b.CalcNormalizedSpectrum(nu_spectrum, uncertainties = uncertainties)

    def CalcDecayChain(self, ZAI: int, tbegin: float = 0., tend: float = np.inf) -> list:
        """
        Calculate fractions of all isotopes in ZAI beta decay chain decaying in specified time window

        :param ZAI: top-of-chain isotope ZAI
        :param begin: begining of the window (in seconds)
        :param end: end of the window (in seconds)
        :return: [(ZAI, f), ...] fractions decaying throughout chain
        """

        res = []
        f0 = 1.
        f1 = 1.

        while ZAI in self.isots:
            HL = self.isots[ZAI].HL
            f0 *= 1. - 2**(-tbegin/HL)
            f1 *= 1. - 2**(-tend/HL)
            res.append((ZAI, f1 - f0))
            ZAI = 10*((ZAI + 1e4)//10)

        return res

    def CalcDecays(self, ZAIs: list, tbegin: float = 0., tend: float = np.inf) -> dict:
        """
        Calculate combined decay chain {ZAI: n, ...} from multiple initial [(ZAI, n), ...]
        """
        res = {}
        for ZAI0, n in ZAIs:
            for ZAI,f in self.CalcDecayChain(ZAI0, tbegin, tend):
                res[ZAI] = res.setdefault(ZAI, 0.) + n*f
        return res

    def SumBranches(self, bset, Es):
        """Sum beta branch spectra with weights {ZAI: w, ...} or [(ZAI, w), ...]"""
        if type(bset) == type({}): bset = bset.items()
        res = Es * 0.
        for ZAI, w in bset:
            for f, df, b in self.isots[ZAI].branches:
                res += w * f * b(Es)
        return res



##############
##############

        # Look for the decay daughters; if they are also beta-unstable, continue to the next generation
        #while self.daughterZAI(generation) in betaSpectraDB.isots.keys():
        #    currentistp = betaSpectraDB.isots[self.daughterZAI(generation)]
        #    HLs.append(currentistp.HL)

        # delta_rate = self.decay_fraction(end, HLs) - self.decay_fraction(begin, HLs) # decay fraction in current window
        # self.decay_chain_spectrum += delta_rate * currentistp.spectrum
        # self.decay_chain_uncertainty += (delta_rate * currentistp.uncertainty)**2

       #     generation += 1

        #self.decay_chain_uncertainty = np.sqrt(self.decay_chain_uncertainty)

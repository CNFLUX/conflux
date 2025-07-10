# Copyright 2025 Lawrence Livermore National Security, LLC. See the top-level NOTICE file for details.
# Author: Xianyi Zhang

# SPDX-License-Identifier: MIT

"""Public modules."""
import numpy as np
import xml.etree.ElementTree as ET
from copy import deepcopy
from tqdm import tqdm
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.linalg import expm

"""CONFLUX modules."""
from conflux.config import CONFLUX_DB
from conflux.Basic import Spectrum, integrate_trapezoid
from conflux.bsg.Constants import ELECTRON_MASS_MEV, NATURAL_LENGTH
from conflux.bsg.SpectralFunctions import (phase_space, 
                                           fermi_function, 
                                           finite_size_L0, 
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

def decay_rhs(t, N, lambdas):
    """Bateman ODEs for a linear chain."""
    dN = np.empty_like(N)
    dN[0] = -lambdas[0] * N[0]
    for i in range(1, N.size):
        dN[i] =  lambdas[i-1] * N[i-1] - lambdas[i] * N[i]
    return dN

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
        #* finite_size_L0_simple(W, **p)
        * recoil_gamow_teller(p["We"], **p)
        * recoil_Coulomb_gamow_teller(p["We"], **p)
        * atomic_screening(p["We"], **p)
    )

    if isNu: result *= radiative_correction_neutrino(p["We"], **p)
    else:    result *= radiative_correction(p["We"], **p)

    if p['L'] == 0:
        result *= shape_factor_gamow_teller(p["We"], **p)
    else:
        result *= shape_factor_unique_forbidden(p["We"], **p)

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

class BetaBranch(Spectrum):
    """Beta decay branch with spectrum"""
    
    Z: int
    """The Atomic number of the mother isotope"""
    A: int
    """The Atomic mass of the mother isotope"""
    I: int
    """The isomeric state this decay mode"""
    Q: float
    """Q value of the decay (MeV)"""
    ZAI: int
    """The unique identity of the isotope (ZAI = Z*1e4+A*10+I)"""
    E0: float
    """End-point energy of this branch, unique for each individual beta-unstable isotope"""
    sigma_E0: float
    """Uncertainty of the end-point energy"""
    frac: float
    """The fraction of this branch to the entire decay."""
    sigma_frac: float
    """The uncertainty of the branching fraction."""
    forbiddenness: int
    """Type of forbiddeness """
    numass: float
    """neutrno mass (MeV), by default 0"""
    mixing: float
    """mixing of nonzero neutrino mass, by default 0"""
    Parameters: dict
    """Beta decay function parameters"""
    corr: dict 
    """Correlation vector of this beta decay branch to other branches of the same isotope"""
    cov: dict 
    """Covariance vector of this beta decay branch with other branches of the same isotope"""
    
    def __init__(self, Z, A, I, Q, E0, sigma_E0, frac, sigma_frac,
                forbiddenness=0, bAc=4.7, xbins = np.arange(0, 20, 0.1),
                numass = 0, mixing = 0, custom_func=None):
        
        Spectrum.__init__(self, xbins)
        
        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q
        self.ZAI = Z*1e4 + A*10 + I

        self.E0 = E0
        self.sigma_E0 = sigma_E0
        self.frac = frac
        self.sigma_frac = sigma_frac

        self.forbiddenness = forbiddenness
        
        self.bAc = bAc
        
        self.numass = numass
        self.mixing = mixing
        
        self.custom_func = custom_func
        #Add the parameters of this branch to a dictionary to be passed onto one of the
        #functions that generates spectra from theory (See neutrino/electron above)
        
        self.Parameters = {
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
        
        self.corr = {E0: 1}  # correlation with other branches of the same isotope
        self.cov = {E0:self.sigma_frac**2} # Set the covariance diagonal element to the square of the branch fraction uncertainty

    # display the vital info of branch
    def Display(self):
        """Display vital branch info."""
        #A little self explanatory as to what's happening
        print("Branch E0 = " + str(self.E0)+ "+\-" + str(self.sigma_E0)
            + ", frac = " + str(self.frac) + "+\-"+str(self.sigma_frac)
            + f", forbiddenness = {self.forbiddenness}")
        
    def UpdateParams(self):
        self.Parameters= {
            'Z': self.Z,
            'A': self.A,
            'W0': self.E0/ELECTRON_MASS_MEV + 1,
            'R': getEltonNuclearRadius(self.A) * 1e-15 / NATURAL_LENGTH,
            'L': self.forbiddenness,
            'c': 1.0,
            'b': self.bAc*self.A,
            'd': 0.0,
            'Lambda': 0.0,
            'l': screening_potential(self.Z)
        }

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
        if self.E0 == otherBranch.E0:
            return
        #Otherwise, set the correlation of the other branch at that endpoint energy to the correlation of this branch
        self.corr[otherBranch.E0] = correlation
        #Also, calculate the covariance of the other branch from the fractoinal uncertainties of the two branches and their correlation
        self.cov[otherBranch.E0] = (self.sigma_frac * correlation
                                    * otherBranch.sigma_frac)

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
        self.Parameters = {
            'Z': self.Z,
            'A': self.A,
            'W0': self.E0/ELECTRON_MASS_MEV + 1,
            'R': getEltonNuclearRadius(self.A) * 1e-15 / NATURAL_LENGTH,
            'L': self.forbiddenness,
            'c': 1.0,
            'b': self.bAc*self.A,
            'd': 0.0,
            'Lambda': 0.0,
            'l': screening_potential(self.Z)
        }
 
        if (nu_spectrum == True):
            function = lambda e: ((1-self.mixing)*neutrino(e, self.Parameters)
                                  +self.mixing*neutrino(e, self.Parameters, numass=self.numass))
        else:
            function = lambda e: ((1-self.mixing)*electron(e, self.Parameters)
                                  +self.mixing*electron(e, self.Parameters, numass=self.numass))
            
        if (self.custom_func!=None):
            function = lambda e: ((1-self.mixing)*self.custom_func(e, self.Parameters)
                                  +self.mixing*self.custom_func(e, self.Parameters, numass=self.numass))

        return function(e)

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
        newParameters = deepcopy(self.Parameters)
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

    def CalcNormalizedSpectrum(self, nu_spectrum=False):
        """
        Calculated unit-normalized beta/neutrino spectrum shape (at `xbins` coordinates)
        
        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to False
        :type nu_spectrum: bool, optional
        """

        self.spectrum = self.BetaSpectrum(self.xbins, nu_spectrum)
        self.uncertainty = np.array([self.SpectEndpointUncertMC(x, nu_spectrum) for x in self.xbins])

        # Normalization spectrum (calculate full range if not already covered)
        if self.xbins[0] <= 0 and self.xbins[-1] >= self.E0:
            norm = integrate_trapezoid(self.xbins, self.spectrum)
        else:
            nrmx = np.arange(0., self.E0, self.E0/50.)
            nrmspec = self.BetaSpectrum(nrmx, nu_spectrum)
            norm = integrate_trapezoid(nrmx, nrmspec)

        if not norm > 0:
            self.spectrum = np.zeros(self.nbin)
        else:
            self.spectrum /= norm
            self.uncertainty /= norm

class BetaIstp(Spectrum):
    """Collection of all beta decay branches from one isotope"""
    
    ZAI : int
    """Top-of-chain isotope (ZAI = Z*1e4+A*10+I)"""
    Z: int
    """The Atomic number of the mother isotope"""
    A: int
    """The Atomic mass of the mother isotope"""
    I: int
    """The isomeric state this decay mode"""
    HL: float
    """Half life of this isotope (s)"""
    Q: float
    """Q value of the decay (MeV)"""
    name: str
    """Name of the isotope"""
    missing: bool = False
    """Uncertainty of the end-point energy"""
    branches: dict
    """Dictionary of decay brances, keys (float) are endpoint energies, values are :class:`conflux.BetaEngine.BetaBranch` """
    numass: float
    """neutrno mass (MeV), by default 0"""
    mixing: float
    """mixing of nonzero neutrino mass, by default 0"""
    MaxBranch: "BetaEngine.BetaBranch"
    """The beta decay branch with the largest fraction"""
    
    
    def __init__(self, Z, A, I, Q, HL, name, xbins=np.arange(0, 20, 0.1), numass=0, mixing=0):
        """Constructor method."""
        Spectrum.__init__(self, xbins)
        
        self.HL = HL # top-of-chain isotope half-life (seconds)

        self.Z = Z
        self.A = A
        self.I = I
        self.ZAI = self.daughterZAI(0)
        self.Q = Q
        self.name = name
        self.missing=False
        self.branches = {}
        
        self.numass = numass
        self.mixing = mixing

    def AddBranch(self, branch):
        """
        Add beta branch with a new endpoint energy to this isotope.
        
        :param branch: the branch to be added
        :type branch: :class: `conflux.BetaEngine.BetaBranch`
        """
        self.branches[branch.E0] = branch


    def EditBranch(self, defaultE0, **kwargs):
        """
        Add or edit branches to the isotope with analyzer's assumptions.
        
        :param defaultE0: Existing end point energy of the missing branch
        :type defaultE0: float
        :param **kwargs: AAdditional keyword arguments.

        Keyword Args:
            sigma_E0 (float): The uncertainty of the end-point energy.
            fraction (float): The fraction that this branch contributes to the total isotopic spectrum.
            sigma_frac (float): The uncertainty of fraction.
            forbiddenness (int): The type of forbidden/allowed transition.
            bAc (float): weak magnetism correcion.
            custom_func (method): customized beta function, needs to be structued similarly as :meth:`neutrino` or :meth:`electron`, defaults to None.
        """
        # if sigma_E0 > E0:
        #     sigma_E0 = E0
        branchexist = defaultE0 in self.branches.keys()

        # setting up default values
        sigma_E0 = self.branches[defaultE0].sigma_E0 if branchexist else 0
        fraction = self.branches[defaultE0].frac if branchexist else 0
        sigma_frac = self.branches[defaultE0].sigma_frac if branchexist else 0
        forbiddenness = self.branches[defaultE0].forbiddenness if branchexist else 0
        bAc = 4.7
        custom_func = self.branches[defaultE0].custom_func if branchexist else None

        for key, value in kwargs.items():
            if key == 'E0': 
                # when E0 is given, replace the default E0
                if branchexist: 
                    del self.branches[defaultE0]
                    branchexist = False
                defaultE0 = value
            if key == 'sigma_E0':
                sigma_E0 = value
            if key == 'fraction':
                fraction = value
            if key == 'forbiddenness':
                forbiddenness = value
            if key == 'bAc':
                bAc = value
            if key == 'custom_func':
                custom_func = custom_func

            if not branchexist:
                # TODO need to remove  the original E0 when a new E0 is given
                self.branches[defaultE0] = BetaBranch(self.Z, 
                                                      self.A, 
                                                      self.I, 
                                                      self.Q, 
                                                      defaultE0, 
                                                      sigma_E0, 
                                                      fraction, 
                                                      sigma_frac,
                                                      forbiddenness,
                                                      bAc=bAc, 
                                                      xbins=self.xbins,
                                                      custom_func=custom_func,
                                                      numass=self.numass,
                                                      mixing=self.mixing)

            elif hasattr(self.branches[defaultE0], key):
                setattr(self.branches[defaultE0], key, value)
                self.branches[defaultE0].UpdateParams()

    def MaxBranch(self):
        """
        Get the beta branch with the maximum end-point energy.
        
        :return: the maximum end-point energy betabranch
        :rtype: :class:`conflux.BetaEngine.BetaBranch`

        """
        Emax = max(self.branches)
        return self.branches[Emax]

    def CalcCovariance(self, GSF=True):
        """
        Calculate the covaraince matrix for all the beta branches of this isotope.
        
        :param GSF: Determine whether to consider the anti-correlation between the Ground-State-decay branching Fraction (GSF) and other branches, defaults to True
        :type GSF: bool, optional

        """
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

    def SumSpectra(self, nu_spectrum=True, branchErange=[0, 20.]):
        """
        Sum the beta/antineutrino spectra of all branches of this isotope.
        
        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to True
        :type nu_spectrum: bool, optional
        :param branchErange: set the lower and upper limit of the endpoint energy of branches to be summed, defaults to [0, 20.]
        :type branchErange: two-element list, optional
        """

        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)
        self.spectUnc = np.zeros(self.nbin) # theoretical uncertainty
        self.branchUnc = np.zeros(self.nbin)
        totalUnc = np.zeros(self.nbin)

        # calculate the total uncertaintty and append spectra
        for E0i, branchi in self.branches.items():
            if(E0i < branchErange[0] or E0i > branchErange[1]):
                continue

            branchi.CalcNormalizedSpectrum(nu_spectrum)
            si = branchi.spectrum
            fi = branchi.frac

            self.spectrum += si * fi
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

    def Display(self):
        """Display vital isotope property and branch information."""
        print('Beta isotope: '+self.name+', ZAI = '+str(self.ZAI)+', Q = '
            +str(self.Q)+" MeV, " +str(len(self.branches))+" branches")
        for E0, branch in self.branches.items():
            branch.Display()

    def decay_time_adjust(self, begin, end):
        """
        Calculate the percentage of isotope decayed in the given time window.
        
        :param begin: the begining of the window (s)
        :type begin: float
        :param end: the end of the window (s)
        :type end: float
        :return: The fraction of isotope decayed (from 0 to 1)
        :rtype: float
        """
        return 2**(-begin/self.HL) - 2**(-end/self.HL)

    def daughterZAI(self, gen):
        """
        ZAI identifier for decay daughter at generation `gen`.
        
        :param gen: daughter generation
        :type gen: int
        :return: The ZAI of the decay daughter
        :rtype: int
        """
        return (self.Z+gen)*1e4 + self.A*10 + (0 if not gen else self.I)
    
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
    
    def CalcDecayChainOld(self, betaSpectraDB, begin, end):
        """
        Calculate the total spectrum of a beta-decay chain in a selected window.
        Assumes beta decays are 100% of allowed decays in chain.
        
        :param betaSpectraDB: the spectrum database that saves all relavant spectra
        :type betaSpectraDB: :class:`conflux.BetaEngine.BetaEngine`
        :param begin: begining of the window (in seconds)
        :type begin: float
        :param end: end of the window (in seconds)
        :type end: float
        :return: summed, decay rate adjusted spectrum and uncertainty in the calculated window
        :rtype: :class:`numpy.array`
        """

        generation = 0
        HLs = []
        self.decay_chain_spectrum = np.zeros(len(betaSpectraDB.xbins))
        self.decay_chain_uncertainty = np.zeros(len(betaSpectraDB.xbins))
        
        # Look for the decay daughters; if they are also beta-unstable, continue to the next generation
        while self.daughterZAI(generation) in betaSpectraDB.istplist.keys():
            currentistp = betaSpectraDB.istplist[self.daughterZAI(generation)]
            HLs.append(currentistp.HL)
            
            delta_rate = self.decay_fraction(end, HLs) - self.decay_fraction(begin, HLs) # decay fraction in current window
            self.decay_chain_spectrum += delta_rate * currentistp.spectrum
            self.decay_chain_uncertainty += (delta_rate * currentistp.uncertainty)**2
            
            generation += 1
            
        self.decay_chain_uncertainty = np.sqrt(self.decay_chain_uncertainty)
        
        return self.decay_chain_spectrum, self.decay_chain_uncertainty
    
    def CalcDecayChain(self, betaSpectraDB, time):
        """
        Calculate the total spectrum of a beta-decay chain in a selected window.
        Assumes beta decays are 100% of allowed decays in chain.
        
        :param betaSpectraDB: the spectrum database that saves all relavant spectra
        :type betaSpectraDB: :class:`conflux.BetaEngine.BetaEngine`
        :param time: the time stamp 
        :type time: float
        :return: summed, decay rate adjusted spectrum and uncertainty in the calculated window
        :rtype: :class:`numpy.array`
        """
        generation = 0
        isotopes = []
        HLs = []

        self.decay_chain_spectrum = np.zeros(len(betaSpectraDB.xbins))
        self.decay_chain_uncertainty = np.zeros(len(betaSpectraDB.xbins))
        
        # Look for the decay daughters; if they are also beta-unstable, continue to the next generation
        while self.daughterZAI(generation) in betaSpectraDB.istplist.keys():
            currentistp = betaSpectraDB.istplist[self.daughterZAI(generation)]
            isotopes.append(self.daughterZAI(generation))
            HLs.append(currentistp.HL)
            
            generation += 1
        
        lambdas = np.log(2)/np.array(HLs)
        k = len(lambdas)
        A = np.zeros((k, k))
        for i, lam in enumerate(lambdas):
            A[i, i] = -lam
            if i > 0:
                A[i, i-1] = lambdas[i-1]

        N0 = np.zeros(k)
        N0[0] = 1
        N_t = expm(A * time).dot(N0)

        decayrates = lambdas*N_t

        for i in range(len(isotopes)):
            ZAI = isotopes[i]
            decayrate = decayrates[i]
            self.decay_chain_spectrum += decayrate * betaSpectraDB.istplist[ZAI].spectrum
            self.decay_chain_uncertainty += (decayrate * betaSpectraDB.istplist[ZAI].uncertainty)**2

        self.decay_chain_uncertainty = np.sqrt(self.decay_chain_uncertainty)

        return self.decay_chain_spectrum, self.decay_chain_uncertainty


# BetaEngine tallys beta branches in the betaDB and calculate theoretical beta 
# spectra of all tallied branches
# if inputlist is not given, load the entire betaDB from the default betaDB
class BetaEngine:
    """
    BetaEngine tallys beta branches in the betaDB and calculate theoretical beta spectra of all tallied branches.
    """
    inputlist: list[int] = None
    """A list of isotopes. If the inputlist is not given, load the entire betaDB from the default betaDB."""
    istplist: dict
    """A dictionary of isotopes. istplist contain keys as the ZAI number of the isotope and values being :class:`conflux.BetaEngine.BetaIstp`"""
    targetDB: str = CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml"
    """The file name of beta decay data base, defaults to CONFLUX_DB+`/betaDB/ENSDFbetaDB2.xml'"""
    xbins: np.ndarray
    """The spectrum range and binning, defaults to np.arange(0, 20, 0.1) (MeV)"""
    custom_func: callable = None
    """Customized beta function, needs to be structued similarly as :meth:`conflux.BetaEngine.neutrino` or :meth:`conflux.BetaEngine.electron`, defaults to None"""
    numass: float = 0
    """To calculate the spectrum with non-zero neutrino mass, give the neutrino mass in the MeV unit."""
    mixing: float = 0
    """To calcualte spectrum with non-zero neutrino mass, provide the mixing (0-1) of the neutrino mass state."""
    
    
    def __init__(self, 
                 inputlist=None, 
                 targetDB=CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml",
                 xbins=np.arange(0, 20, 0.1),
                 custom_func=None,
                 numass=0,
                 mixing=0):
        """Constructor method."""
        self.inputlist = inputlist
        self.istplist = {}
        self.xbins = xbins
        self.custom_func=custom_func
        
        self.numass = numass
        self.mixing = mixing

        self.LoadBetaDB(targetDB)   # loadBetaDB automatically
        
    def LoadBetaDB(self, targetDB=CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml", missingBranch = 3):
        """
        Load default or input betaDB to obtain beta decay informtion. A customed DB must follow the same format as the default DB.
        
        :param targetDB: The file name of beta decay data base, defaults to CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml"
        :type targetDB: str, optional
        :param missingBranch: Determine how many branches there are in a missing isotope, defaults to 3
        :type missingBranch: int, optional

        """
        # empty the existing isotope list before filling loading the DB
        self.istplist= {}
        
        useInputList = True # test if the engine is defined with an inputlist
        if self.inputlist == None:
            print("Loading all beta data from the default betaDB...")
            self.inputlist = []
            useInputList = False

        print("Searching DB: "+targetDB+"...")
        print("Loading spectra of beta branches...")

        tree = ET.parse(targetDB)
        root = tree.getroot()
        for isotope in root:
            ZAI = int(isotope.attrib['isotope'])
            if ZAI == 10: continue # ignoring single neutrons
            Q = float(isotope.attrib['Q'])
            HL = float(isotope.attrib['HL'])
            name = isotope.attrib['name']

            # if input list is not given, include all isotopes
            if not useInputList:
                self.inputlist.append(ZAI)
            

            if (ZAI in self.inputlist):
                Z = int(ZAI/10000)
                A = int(ZAI%10000/10)
                I = int(ZAI%10)

                # if same ZAI reappear it is tagged as an isomer (unmarked in
                # the original DB)
                if ZAI in self.istplist:
                    I += 1
                    ZAI += 1

                betaIstp = BetaIstp(Z, A, I, Q, HL, name, xbins=self.xbins)

                # Adding missing branches below
                if len(isotope) < 1:
                    betaIstp.missing = True
                    for i in range(missingBranch):
                        misfrac = 1./missingBranch
                        misE0 = betaIstp.Q/missingBranch*(i+1)
                        betaIstp.EditBranch(misE0, E0=misE0, fraction=misfrac)
                    self.istplist[ZAI] = betaIstp
                    continue

                # some isotopes contain summed branch fraction greater than 1
                # prepare to normalize
                fracsum = 0
                for branch in isotope:
                    fraction = float(branch.attrib['fraction'])
                    fracsum += fraction

                # actually assign values from database to branches
                for branch in isotope:
                    E0 = float(branch.attrib['end_point_E'])
                    if E0<=0:
                        continue
                    sigma_E0 = float(branch.attrib['sigma_E0'])
                    spin_par_changes = [n.strip() for n in branch.attrib['dJpi'].split(',')]
                    # spin_par_changes = str(branch.attrib['dJpi'])

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
                    if fracsum > 1:
                        fraction /= fracsum

                    betaBranch = BetaBranch(Z, A, I, Q, E0, sigma_E0, fraction,
                                            sigma_frac, forbiddenness,
                                            xbins=self.xbins,
                                            custom_func=self.custom_func,
                                            numass=self.numass,
                                            mixing=self.mixing)
                    betaIstp.AddBranch(betaBranch)

                if betaIstp.branches:
                    self.istplist[ZAI] = betaIstp

    def GetBetaIstp(self, Z, A, I=0):
        """
        Get the BetaIstp object by giving the ZAI number
        
        :param Z: Proton number of the isotope
        :type Z: int
        :param A: Atomic mass number of the isotope
        :type A: int
        :param I: the Ith isomeric state, defaults to 0
        :type I: int, optional
        :return: the specified beta isotope object
        :rtype: :class:`conflux.BetaEngine.BetaIstp`

        """
        ZAI = Z*1e4 + A*10 + I
        return self.istplist[ZAI]
        
    def CalcBetaSpectra(self, nu_spectrum=True,
                        branchErange=[0.0, 20.0], GSF=False, silent=False):
        """
        Calculates beta or neutrino spectra of listed/default beta decaying isotopes.
        
        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to True
        :type nu_spectrum: bool, optional
        :param branchErange: set the lower and upper limit of the endpoint energy of branches to be summed, defaults to [0, 20.]
        :type branchErange: two-element list, optional
        :param GSF: Determine whether to consider the anti-correlation between the Ground-State-decay branching Fraction (GSF) and other branches, defaults to True
        :type GSF: bool, optional
        :param silent: whether to disable the tqdm output, defaults to False
        :type silent: bool, optional

        """
    
        for ZAI,betaIstp in tqdm(self.istplist.items(),
                                 desc=f"Calculating beta/neutrino spectra of {len(self.istplist)} isotopes",
                                 disable=silent):
            if betaIstp.Q < branchErange[0] or betaIstp.Q > branchErange[1]:
                continue
            betaIstp.CalcCovariance(GSF=GSF)
            betaIstp.SumSpectra(nu_spectrum, branchErange)
            
    def SaveToFile(self, filename):
        spect_data = {}
        spect_data["xbins"] = self.xbins
        for ZAI,betaIstp in self.istplist.items():
            spect_data[f"{ZAI}"] = betaIstp.spectrum
            spect_data[f"{ZAI}_unc"] = betaIstp.uncertainty
            
        df = pd.DataFrame(spect_data)
        df.to_csv(filename, index=False)
            
    def LoadFile(self, filename):
        df = pd.read_csv(filename)
        self.xbins = df["xbins"].values
        
        for col in df.columns:
            if col == "xbins":
                continue
            elif not col.endswith("_unc"):
                key = int(col)  # remove "_unc"
                Z = key/1e4
                A = (key-Z*1e4)/10
                I = key-Z*1e4-A*10
                # newistp = BetaIstp(Z, A, I, Q=0, HL=0, name="", xbins=self.xbins)
                self.istplist[key].spectrum=df[col].values
            if col.endswith("_unc"):
                key = int(col[:-4])
                self.istplist[key].unsertinaty=df[col].values        
                # self.istplist[key]=newistp
                
        print(f"Loaded spectra and uncertainties from {filename}")
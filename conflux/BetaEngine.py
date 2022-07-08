import sys
import os
import argparse
import numpy as np
from scipy import special, interpolate, integrate
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from copy import deepcopy
import pkg_resources
import timeit

from conflux.bsg.Constants import *
from conflux.bsg.SpectralFunctions import *
from conflux.bsg.Functions import *
from conflux.bsg.Screening import *

#########################################
# Final neutrino and antineutrino spectra

def electron(ebeta, p):
    result = 0.
    W = ebeta/ELECTRON_MASS_MEV + 1
    result = (phase_space(W, **p)
            *fermi_function(W, **p)
            *finite_size_L0(W, **p)
            *recoil_gamow_teller(W, **p)
            *radiative_correction(W, **p)
            *recoil_Coulomb_gamow_teller(W, **p)
            *atomic_screening(W, **p)
            )

    if p['L'] == 0:
        result *= shape_factor_gamow_teller(W, **p)
    else:
        result *= shape_factor_unique_forbidden(W, **p)

    return result

def neutrino(enu, p):
    result = 0.
    W0 = p['W0']
    Wv = W0-enu/(ELECTRON_MASS_MEV*1.0) #enu/ELECTRON_MASS_MEV + 1
    result = (phase_space(Wv, **p)
            *fermi_function(Wv, **p)
            *finite_size_L0(Wv, **p)
            *recoil_gamow_teller(Wv, **p)
            *radiative_correction_neutrino(Wv, **p)
            *recoil_Coulomb_gamow_teller(Wv, **p)
            *atomic_screening(Wv, **p)
            )

    if p['L'] == 0:
        result *= shape_factor_gamow_teller(Wv, **p)
    else:
        result *= shape_factor_unique_forbidden(Wv, **p)

    return result

def integral(nu_spectrum, p, x_low, x_high):
    erg = 0.0
    xh = 0.0
    x_low = max([x_low, 1e-6])

    if (x_low >= p.e0-1e-6):
        return 0
    if (x_high > p.e0):
        xh = p.e0
    else:
        xh = x_high

    if (nu_spectrum == True):
        function = lambda x: neutrino(x, p)
    else:
        function = lambda x: electron(x, p)

    #function(x_low)
    if (x_low <= p.thresh and xh >= p.thresh): # wierd error because of integrater when difference is smaller than a value

        ergl = [0.0, 0.0]
        ergh = [0.0, 0.0]

        # a quick integration mathod through trapezoid algorithm
        mid = (x_low+p.thresh)/2
        ergl = abs(p.thresh-x_low)*function(mid)
        mid = (xh+p.thresh)/2
        ergh = abs(p.thresh-xh)*function(mid)
        erg = ergl+ergh

        return erg

    mid = (xh+x_low)/2
    erg = abs(xh-x_low)*function(mid)

    return erg

    #print(x_low, p.thresh, x_high, erg[0])

# BetaBranch class to save the isotopic information
class BetaBranch:
    def __init__(self, Z, A, I, Q, E0, sigma_E0, frac, sigma_frac, forbiddeness=0, bAc=4.7):
        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q
        self.ZAI=Z*1e4+A*10+I

        self.E0 = E0
        self.sigma_E0 = sigma_E0
        self.frac = frac
        self.sigma_frac = sigma_frac

        self.forbiddeness = forbiddeness

        self.Parameters = {
            'Z': Z,
            'A': A,
            'W0': E0/ELECTRON_MASS_MEV + 1,
            'R': getEltonNuclearRadius(A) * 1e-15 / NATURAL_LENGTH,
            'L': forbiddeness,
            'c': 1.0,
            'b': bAc*A,
            'd': 0.0,
            'Lambda': 0.0,
            'l': screening_potential(Z)
        }

        #self.Parameters = Parameters_t(Z = self.Z+1, A = self.A, e0=self.E0, WM=self.WM, ftype=self.forbiddeness)

        self.corr = {E0:1}  # correlation with other branches of the same isotope
        self.cov = {E0:self.sigma_frac**2}

    # display the vital info of branch
    def Display(self):
        print("Branch E0 = " +str(self.E0)+"+\-"+str(self.sigma_E0)+", frac = "+str(self.frac)+"+\-"+str(self.sigma_frac))

    # set correlation of this branch fraction and all another branches
    def SetCovariance(self, otherBranch, correlation):
        if self.E0 == otherBranch.E0:
            return
        self.corr[otherBranch.E0] = correlation
        self.cov[otherBranch.E0] = self.sigma_frac*correlation*otherBranch.sigma_frac

    # beta spectrum shape as function of energy
    def BetaSpectrum(self, x, nu_spectrum=False):
        Parameters = deepcopy(self.Parameters)
        rangeCorrect = x <= self.E0 # prevent out-of-range variable to output insane results

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, Parameters)
        else:
            function = lambda x: electron(x, Parameters)

        result = function(x)
        result = np.nan_to_num(result, nan=0.0)
        return result*rangeCorrect

    # calculate the spectrum uncertainty
    def SpectUncert(self, x, nu_spectrum=False):

        numbers = 5
        E0range = np.linspace(self.E0-self.sigma_E0, self.E0+self.sigma_E0, numbers)
        newParameters = deepcopy(self.Parameters)
        newParameters['W0'] = E0range/ELECTRON_MASS_MEV+1

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, newParameters)
        else:
            function = lambda x: electron(x, newParameters)

        fE0 = function(x)
        fE0 = np.nan_to_num(fE0, nan=0.0)
        # Commenting out unknown function
        # if (fE0.all() < 1e-8):
        #     return 0

        grad_E0 = np.gradient(fE0)

        return grad_E0[int(numbers/2)]*self.sigma_E0

    # calculate the spectrum uncertainty with MC sampling
    def SpectUncertMC(self, x, nu_spectrum=False, samples = 50):

        E0range = np.random.normal(self.E0, self.sigma_E0, samples)
        newParameters = deepcopy(self.Parameters)
        newParameters['W0'] = E0range/ELECTRON_MASS_MEV+1

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, newParameters)
        else:
            function = lambda x: electron(x, newParameters)

        fE0 = function(x)
        fE0 = np.nan_to_num(fE0, nan=0.0)
        # Commenting out unknown function
        # if (fE0.all() < 1e-8):
        #     return 0
        return np.std(fE0)

    # bined beta spectrum
    def BinnedSpectrum(self, nu_spectrum=False, binwidths=0.1, spectRange=[-1.0, 20.0]):
        bins = int(spectRange[1]/binwidths)
        self.result = np.zeros(bins)
        self.uncertainty = np.zeros(bins)

        lower = spectRange[0]
        if (lower > self.E0):
            return 1
        if (lower<0):
            lower=binwidths/2.0

        # integrating each bin
        # startTiming = timeit.default_timer()

        for k in range(0, bins):
            x_low = lower
            x_high = lower+binwidths
            if x_high > self.E0:
                x_high = self.E0
            self.result[k] = abs(x_high-x_low)*self.BetaSpectrum((x_low+x_high)/2, nu_spectrum)
            #self.uncertainty[k] = abs(x_high-x_low)*(self.SpectUncert(x_low, nu_spectrum)+self.SpectUncert(x_high, nu_spectrum))/2
            #gradUnc = self.uncertainty[k]
            #print("gradUnc", gradUnc)
            self.uncertainty[k] = abs(x_high-x_low)*self.SpectUncertMC((x_low+x_high)/2, nu_spectrum)
            if x_high == self.E0:
                break
            #MCUnc = self.uncertainty[k]
            #print("MCUnc", MCUnc)
            lower+=binwidths

        # endTiming = timeit.default_timer()
        # runTime = endTiming-startTiming
        # print("runtime", runTime)

        # normalizing the spectrum
        norm = self.result.sum()
        #print(self.result, self.uncertainty, norm*binwidths)
        if norm <=0:
            self.result =np.zeros(bins)
        else:
            self.result /= norm*binwidths
            self.uncertainty /= norm*binwidths
            #print(self.result, self.uncertainty, norm*binwidths)
        return 0

# class to save isotope information, including Z A I Q and beta brances
class BetaIstp:
    def __init__(self, Z, A, I, Q, name):
        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q
        self.name = name
        self.ZAI=Z*1e4+A*10+I
        self.branch={}
        self.missing=False

    def AddBranch(self, branch):
        """
        Add beta branch to this isotope
        Returns:
            None
        """
        self.branch[branch.E0] = branch

    def EditBranch(self, E0, fraction, sigma_E0 = 0., sigma_frac = 0., forbiddeness = 0):
        """
        Add or edit branches to the isotope with analyzer's assumptions
        Parameters:
            E0 (float): assumed end point energy of the missing branch
            fraction (float): assumed fraction of the edited branch
            sigma_E0 (float): 1-sigma error of the input E0, default as 0
            sigma_frac (float): 1-sigma error of the input franction, default as 0
        Returns:
            None
        """
        self.branch[E0] = BetaBranch(self.Z, self.A, self.I, self.Q, E0, sigma_E0, fraction, sigma_frac, forbiddeness)

    def MaxBranch(self):
        """
        Get the maximum E0 branch of the isotope
        Returns:
            self.branch(Emax) (BetaBranch)
        """
        Emax = max(self.branch)
        return self.branch[Emax]

    def CalcCovariance(self, GSF=True):
        """
        Calculate the covaraince matrix for all the beta branches of this isotope
        Parameters:
            GSF (boolean): determine whether to calculate covaraince matrix with ground state feeding
        Returns:
            None
        """
        MaxBranch = self.MaxBranch()
        # obtain the ground state branch info
        totalFrac = 0
        GSFrac = MaxBranch.frac
        for E0, branch in self.branch.items():
            totalFrac += branch.frac
        restFrac = totalFrac-GSFrac

        if GSF == False or MaxBranch.E0 != self.Q or restFrac == 0:
            for i, branchi in self.branch.items():
                for j, branchj in self.branch.items():
                    branchi.SetCovariance(branchj, 0)
            return

        # calculate the covariance matrix with ground state anticorrelation
        for i, branchi in self.branch.items():
            for j, branchj in self.branch.items():
                if i == self.Q:
                    correlation = -1*branchj.frac/restFrac
                    branchi.SetCovariance(branchj, correlation)
                    branchj.SetCovariance(branchi, correlation)
                elif j != self.Q:
                    branchi.SetCovariance(branchj, 0)

    def CalcBetaSpectrum(self, nu_spectrum=True, binwidths=0.1, spectRange=[-1.0, 20.0]):
        """
        Calculate the cumulative beta/antineutrino spectrum of all branches
        Parameters:
            nu_spectrum (boolean): Determines if you are calculating a beta spectra or a neutrino spectra
            binwidths (float): The width of the bins used in creating your spectra
            spectRange (list): define the upper and power bounds of the spectrum
        Returns:
            None
        """
        spectLow = spectRange[0] if spectRange[0] > 0 else 0
        spectHigh = spectRange[1]
        bins = int((spectHigh-spectLow)/binwidths)
        self.spectrum=np.zeros(bins)
        self.spectUnc=np.zeros(bins)
        self.branchUnc=np.zeros(bins)
        self.totalUnc=np.zeros(bins)

        for E0,branch in self.branch.items():
            branch.BinnedSpectrum(nu_spectrum, binwidths, spectRange)
            self.spectrum += branch.result*branch.frac
            self.spectUnc += branch.uncertainty*branch.frac

        for E0i, branchi in self.branch.items():
            si = branchi.result
            di = branchi.sigma_frac
            for E0j, branchj in self.branch.items():
                sj = branchj.result
                dj = branchj.sigma_frac
                cov_bij = branchi.cov[E0j]
                sigmab_ij = si*cov_bij*sj
                self.branchUnc += sigmab_ij
                if (E0i==E0j):
                    self.totalUnc += (branchi.uncertainty*branchi.frac)**2 + sigmab_ij
                else:
                    self.totalUnc += sigmab_ij
        self.branchUnc = np.sqrt(self.branchUnc)
        self.totalUnc = np.sqrt(self.totalUnc)

    def Display(self):
        """
        Display isotope property and branch information
        """
        print('Beta isotope: '+self.name+', ZAI = '+str(self.ZAI)+', Q = '+str(self.Q)+" MeV, " +str(len(self.branch))+" branches")
        for E0, branch in self.branch.items():
            branch.Display()

# BetaEngine tallys beta branches in the betaDB and calculate theoretical beta spectra
# of all tallied branches
# if inputlist is not given, load the entire betaDB from the default betaDB
class BetaEngine:
    def __init__(self, inputlist=None, targetDB=None):
        self.inputlist = inputlist
        self.defaultDB = os.environ["CONFLUX_DB"]+"/betaDB/ENSDFbetaDB.xml"

        self.LoadBetaDB(targetDB)   # loadBetaDB automatically

    def LoadBetaDB(self, targetDB=None):
        useInputList = True         # test if the engine is defined with an inputlist
        if self.inputlist == None:
            print("Loading all beta data from the default betaDB...")
            self.inputlist = []
            useInputList = False

        if targetDB == None:
            targetDB = self.defaultDB
        print("Searching DB: "+targetDB+"...")
        print("Loading spectra of beta branches...")

        self.istplist = {}

        tree = ET.parse(targetDB)
        root = tree.getroot()
        for isotope in root:
            ZAI = int(isotope.attrib['isotope'])
            Q = float(isotope.attrib['Q'])
            name = isotope.attrib['name']

            if not useInputList:
                self.inputlist.append(ZAI)

            if (ZAI in self.inputlist):
                #print(str(ZA)+"...")
                Z = int(ZAI/10000)
                A = int(ZAI%10000/10)
                I = int(ZAI%10)

                betaIstp = BetaIstp(Z, A, I, Q, name)

                # Adding missing branches below
                if len(isotope) < 1:
                    betaIstp.missing = True
                    betaIstp.EditBranch(betaIstp.Q, 1)
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
                    sigma_E0 = float(branch.attrib['sigma_E0'])
                    spin_par_changes = (branch.attrib['dJpi'])

                    # converting spin change to forbidden types
                    ftypes = [['0', '1'], ['0-', '1-', '2-'], ['2', '3'], ['3-', '4-'], ['4', '5']]
                    firstftypes = [-10, -11, 10]
                    forbiddeness = 1e3
                    for i in range(len(ftypes)):
                        for j in range(len(spin_par_changes)):
                            if spin_par_changes[j] in ftypes[i] and i < abs(forbiddeness):
                                forbiddeness = -i
                                if i == 1:
                                    forbideness = firstftypes[j]
                                elif spin_par_changes[j] == ftypes[i][-1]:
                                    forbideness = i

                    # assign fraction values to branches
                    # normalize if greater than one
                    fraction = float(branch.attrib['fraction'])
                    sigma_frac = float(branch.attrib['sigma_frac'])
                    if fracsum > 1:
                        fraction /= fracsum
                        sigma_frac /= fracsum

                    betaBranch = BetaBranch(Z, A, I, Q, E0, sigma_E0, fraction, sigma_frac, forbiddeness)
                    betaIstp.AddBranch(betaBranch)
                self.istplist[ZAI] = betaIstp


    def CalcBetaSpectra(self, targetDB = None, nu_spectrum=True, binwidths=0.1, spectRange=[-1.0, 20.0], branchErange=[-1.0, 20.0]):
        """
        Calculates beta spectra of the list of beta-decaying isotopes

        Parameters:
            targetDB (String): the path to the betaSpectra database. if none, use the default database.
            nu_spectrum (boolean): Determines if you are calculating a beta spectra or a neutrino spectra
            binwidths (float): The width of the bins used in creating your spectra.
            spectRange (list): define the upper and power bounds of the spectrum
            branchErange (list): defind the range of interested Q values
        Returns:
            None

        """
        spectLow = spectRange[0] if spectRange[0] > 0 else 0
        spectHigh = spectRange[1]
        self.bins = np.arange(spectLow, spectHigh, binwidths)

        startTiming = timeit.default_timer()
        istpCount = 0
        for ZAI in self.istplist:
            betaIstp = self.istplist[ZAI]
            if betaIstp.Q < branchErange[0] or betaIstp.Q > branchErange[1]:
                continue

            betaIstp.CalcCovariance()
            betaIstp.CalcBetaSpectrum(nu_spectrum, binwidths, spectRange)
            istpCount += 1

        endTiming = timeit.default_timer()
        nBranch = istpCount
        runTime = endTiming-startTiming
        print("Finished calculating beta spectra of "+ str(nBranch) + " isotopes.")
        print("Processing time: "+str(runTime)+" seconds")

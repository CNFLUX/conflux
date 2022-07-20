# fission product and spectra summation engine

# universal modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

# local modules
from conflux.Basic import *
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp

class SumEngine(Spectrum):
    """
    A Class to carry out the Summation (Ab-initio) Mode

    ...

    Attributes
    ----------

    FPYlist : dictionary
        A dictionary of the fission product yield (FPY) of each isotope
    betaSpectraList : dictionary
        A dictionary of the Beta Spectra of each isotope
    betaUncertainty : int list
        A list of the beta Uncertainty's of each isotope
    neutrino : boolean
        A boolean of whether we are looking at the neutrino spectrum or beta spectrum

    Methods
    -------

    Clear():
        Clears all dictionaries associated with the SumEngine

    AddModel(fissionModel, W=1.0):
        method to add fission/non-fissile/non-equilibrium isotopes into the engine

    Normalize():
        Normalizes the FPY and the fission product yield Error

    CalcReactorSpectrum(betaSpectraDB, binwidths=0.1, lower=-1.0, thresh=0.0, erange=20.0):
        Calculates the neutrino Spectrum from the given database with the given energy bins and bounds

    Draw(figname, summing=True, logy=True, frac=False):
        Draws the Reactor Spectrum and saves it as a .png file.

    """
    
    def __init__(self, neutrino=True, binwidths=0.1, spectRange=[-1.0, 20.0]):
        self.FPYlist = {}
        self.betaSpectraList = {}
        self.betaUncertainty = {}
        self.neutrino = neutrino # neutrino or electon spectrum
        
        self.binwidths=binwidths
        self.spectRange=spectRange
        self.bins = int(spectRange[1]/binwidths)
        self.xbins = np.arange(*spectRange, binwidths)
        self.spectrum = np.zeros(self.bins)
        self.uncertainty = np.zeros(self.bins)

    def Clear(self):
        """
            Clears out all associated dictionaries inside the SumEngine

            Parameters:
                None
            Returns:
                None
        """
        self.FPYlist = {}
        self.betaUncertainty = {}
        self.betaSpectraList = {}

    # method to add fission/non-fissile/non-equilibrium isotopes into the engine
    def AddModel(self, fissionModel, W=1.0):
        """
            Adds a reactor model to the summation engine.

            Parameters:
                fissionModel (FissionModel): A Fission Model containing the fission/non-fissile/non-equilibrium isotopes
                W (float) : The weight applied to the reactor model in a multi-reactor engine. The default weight is set to 1.
            Returns:
                None
        """
        for FPZAI in fissionModel.FPYlist:
            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = fissionModel.FPYlist[FPZAI]
                self.FPYlist[FPZAI].y *= W
                self.FPYlist[FPZAI].yerr *= W
            else:
                self.FPYlist[FPZAI].y + fissionModel.FPYlist[FPZAI].y*W
                self.FPYlist[FPZAI].yerr + fissionModel.FPYlist[FPZAI].yerr*W
                self.FPYlist[FPZAI].AddCovariance(fissionModel.FPYlist[FPZAI])

    def NormalizeFP(self):
        """
            Normalizes the fission product yield and the Fission Percent error.

            Parameters:
                None
            Returns:
                None

        """
        self.sum = 0
        for FPZAI in self.FPYlist:
            self.sum += self.FPYlist[FPZAI].y
        for FPZAI in self.FPYlist:
            self.FPYlist[FPZAI].y /= self.sum
            self.FPYlist[FPZAI].yerr /= self.sum

    def CalcReactorSpectrum(self, betaSpectraDB, branchErange=[-1.0, 20.0], processMissing=False):
        """
            Calculates the reactor spectrum based off the fission yield database as well as
            the betaSpectra database.
            Parameters:
                betaSpectraDB (BetaEngine): a dictionary of the spectral information for each beta branch
                binwidths (int): the width of the bins used in running the calculation
                erange (int): upper limit of energy you want to run the reactor for
            Returns:
                None

        """
        
        print("Summing beta spectra...")
        self.modelUnc = np.zeros(self.bins)
        self.yieldUnc = np.zeros(self.bins)
        self.missingBranch = []
        self.missingCount = 0.0
        self.totalYield = 0.0

        for FPZAI in self.FPYlist:
            self.totalYield += self.FPYlist[FPZAI].y
            if FPZAI in betaSpectraDB.istplist:
                
                # process missing branches with assumptions or not
                if not processMissing and betaSpectraDB.istplist[FPZAI].missing:
                    self.missingCount += self.FPYlist[FPZAI].y
                    self.missingBranch.append(FPZAI)
                    continue
                
                self.betaSpectraList[FPZAI] = betaSpectraDB.istplist[FPZAI].spectrum
                self.betaUncertainty[FPZAI] = betaSpectraDB.istplist[FPZAI].totalUnc
                self.spectrum += self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y
                self.yieldUnc += self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].yerr
                self.modelUnc += self.betaUncertainty[FPZAI]*self.FPYlist[FPZAI].y
                self.spectrumUnc = self.yieldUnc+self.modelUnc

            # save the list of missing branches
            else:
                self.missingCount += self.FPYlist[FPZAI].y
                self.missingBranch.append(FPZAI)

        # Uncertainty calculation
        for i in self.FPYlist:
            for j in self.FPYlist:
                if i in self.betaSpectraList and j in self.betaSpectraList:
                    yi = self.FPYlist[i].y
                    yerri = self.FPYlist[i].yerr
                    fi = self.betaSpectraList[i]
                    ferri = self.betaUncertainty[i]

                    yj = self.FPYlist[j].y
                    yerrj = self.FPYlist[j].yerr
                    fj = self.betaSpectraList[j]
                    ferrj = self.betaUncertainty[j]
                    
                    # if covariance matrix were not loaded, make cov diagonal variance
                    if not self.FPYlist[i].cov[j]:
                        sigmay_ij = yerri*yerri if i == j else 0
                    else:
                        sigmay_ij = self.FPYlist[i].cov[j]
                    
                    cov_ij = fi*sigmay_ij*fj
                    
                    if (i==j):
                        self.spectrumUnc += (ferri*yi)**2 + cov_ij*fi**2
                    else:
                        self.spectrumUnc += cov_ij

        self.uncertainty = np.sqrt(self.spectrumUnc)

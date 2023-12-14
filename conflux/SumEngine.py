# fission product and spectra summation engine

# universal modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
from copy import deepcopy

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

    def __init__(self, neutrino=True, xbins=np.arange(0, 20, 0.1)):
        self.FPYlist = {}
        self.betaSpectraList = {}
        self.betaUncertainty = {}
        self.neutrino = neutrino # neutrino or electon spectrum

        self.xbins = xbins
        self.nbin = len(xbins)
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

    # def __del__(self):
    #     print("Summation model cleared")

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
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

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
        inputFPY = deepcopy(fissionModel)
        for FPZAI in inputFPY.FPYlist:
            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = inputFPY.FPYlist[FPZAI]
                self.FPYlist[FPZAI].y *= W
                self.FPYlist[FPZAI].yerr *= W
            else:
                self.FPYlist[FPZAI].y += inputFPY.FPYlist[FPZAI].y*W
                self.FPYlist[FPZAI].yerr += inputFPY.FPYlist[FPZAI].yerr*W
                self.FPYlist[FPZAI].AddCovariance(inputFPY.FPYlist[FPZAI])

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

    def CalcReactorSpectrum(self, betaSpectraDB, branchErange=[-1.0, 20.0], processMissing=False, ifp_begin = 0, ifp_end = 0, modelunc = True):
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
        self.modelUnc = np.zeros(self.nbin)
        self.yieldUnc = np.zeros(self.nbin)
        self.missingBranch = []
        self.missingCount = 0.0
        self.totalYield = 0.0

        for FPZAI in self.FPYlist:
            # get the yield of the fission products
            thisyield = self.FPYlist[FPZAI].y
            yielderr = self.FPYlist[FPZAI].yerr

            self.totalYield += thisyield

            if FPZAI in betaSpectraDB.istplist:

                # for IFP calculation, adjust the decay rate of the target isotope
                # by the rate of isotope that are decayed in the time window
                if (ifp_begin < ifp_end):
                    adjustment = betaSpectraDB.istplist[FPZAI].decay_time_adjust(ifp_begin, ifp_end)
                    thisyield *= adjustment
                    yielderr *= adjustment

                # process missing branches with assumptions or not
                if not processMissing and betaSpectraDB.istplist[FPZAI].missing:
                    self.missingCount += thisyield
                    self.missingBranch.append(FPZAI)
                    continue

                self.betaSpectraList[FPZAI] = betaSpectraDB.istplist[FPZAI].spectrum
                self.betaUncertainty[FPZAI] = betaSpectraDB.istplist[FPZAI].uncertainty
                self.spectrum += self.betaSpectraList[FPZAI]*thisyield

                '''
                obsolete code
                '''
                self.yieldUnc += self.betaSpectraList[FPZAI]*yielderr
                self.modelUnc += self.betaUncertainty[FPZAI]*thisyield
                self.spectrumUnc = self.yieldUnc+self.modelUnc

            # save the list of missing branches
            else:
                self.missingCount += thisyield
                self.missingBranch.append(FPZAI)

        # Uncertainty calculation
        for i in self.FPYlist:
            for j in self.FPYlist:
                if i in self.betaSpectraList and j in self.betaSpectraList:
                    adjustmenti = 1
                    adjustmentj = 1
                    if (ifp_begin < ifp_end):
                        adjustmenti = betaSpectraDB.istplist[i].decay_time_adjust(ifp_begin, ifp_end)
                        adjustmentj = betaSpectraDB.istplist[j].decay_time_adjust(ifp_begin, ifp_end)

                    yi = self.FPYlist[i].y*adjustmenti
                    yerri = self.FPYlist[i].yerr*adjustmenti
                    fi = self.betaSpectraList[i]*adjustmenti
                    ferri = self.betaUncertainty[i]*adjustmenti

                    yj = self.FPYlist[j].y*adjustmentj
                    # yerrj = self.FPYlist[j].yerr*adjustmentj
                    fj = self.betaSpectraList[j]*adjustmentj
                    ferrj = self.betaUncertainty[j]*adjustmentj

                    # if covariance matrix were not loaded, make cov diagonal variance
                    cov_ij = 0
                    if j not in self.FPYlist[i].cov:
                        cov_ij = yerri*yerri if i == j else 0
                    else:
                        cov_ij = self.FPYlist[i].cov[j]

                    sigmay_ij = fi*cov_ij*fj

                    self.uncertainty += sigmay_ij
                    # if (i==j):
                    #     self.uncertainty += sigmay_ij + (ferri*yi)**2
                    # else:
                    #     self.uncertainty += sigmay_ij

        # if allowed, add beta model uncertainty to the result
        if modelunc:
            self.uncertainty += self.modelUnc

        self.uncertainty = np.sqrt(self.uncertainty)

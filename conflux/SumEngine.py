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
        self.FPYlist = {} #Dictionary of Fissile Products
        self.betaSpectraList = {} #Dictionary of beta Spectra for each Fission Product
        self.betaUncertainty = {} #Dictionary of beta Uncertainties for each Fission Product
        self.neutrino = neutrino # neutrino or electon spectrum

        self.xbins = xbins #Energy (in MeV)
        self.nbin = len(xbins) #How fine or course your binning is
        self.spectrum = np.zeros(self.nbin) #The Calculated total antineutrino/beta spectrum
        self.uncertainty = np.zeros(self.nbin) #The Calculated total uncertainties

    #Self explanatory, clears the various dictionaries associated
    #With the Summation Engine.  
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

        #Create copy of the fissionmodel
        inputFPY = deepcopy(fissionModel)
        for FPZAI in inputFPY.FPYlist:
            #Check to see if the fission products in the fission model exists in the Reactor model. 
            
            #If they aren't in the Reactor model, add them to the reactor model, and adjust their
            #Yield and yield error by the weight of the fission mdel.

            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = inputFPY.FPYlist[FPZAI]
                self.FPYlist[FPZAI].y *= W
                self.FPYlist[FPZAI].yerr *= W

            
            #If the fission product already exists in the reactor model, take the yield and yield error
            #that this isotope has in this fission model, and add it to the yield and yield error in the 
            #Reactor model.  
            #Also, Add the covariance of the this fission product from the inputted fission model
            #To the fission product in the Reactor model.  
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
        #Add the yields of all Fission Products
        for FPZAI in self.FPYlist:
            self.sum += self.FPYlist[FPZAI].y
        for FPZAI in self.FPYlist:
            #Divide the yields and errors of each fission product with the total fission yield. 
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
        #initialize total spectrum & uncertainty
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

        
        #Initialize model and Yield Uncertainties, a list of missing branches to the total
        #Contribution, the missing yield, and the total yield. 
        self.modelUnc = np.zeros(self.nbin)
        self.yieldUnc = np.zeros(self.nbin)
        self.missingBranch = []
        self.missingCount = 0.0
        self.totalYield = 0.0

        #Iterate through every fission product in the Reactor Model.        
        for FPZAI in self.FPYlist:
            # get the yield of the fission products
            thisyield = self.FPYlist[FPZAI].y
            yielderr = self.FPYlist[FPZAI].yerr

            #Add the yield of each fission product to the total yield.
            self.totalYield += thisyield

            #If the fission product in the reactor model is also in the Beta Spectra model, do some calculations            
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

                #Pull the beta spectrum and the beta uncertainties, add the product of the Beta Spectra and yield
                #To the total spectrum
                self.betaSpectraList[FPZAI] = betaSpectraDB.istplist[FPZAI].spectrum
                self.betaUncertainty[FPZAI] = betaSpectraDB.istplist[FPZAI].uncertainty
                self.spectrum += self.betaSpectraList[FPZAI]*thisyield

                '''
                obsolete code
                '''
                self.yieldUnc += self.betaSpectraList[FPZAI]*yielderr
                self.modelUnc += self.betaUncertainty[FPZAI]*thisyield
                self.spectrumUnc = self.yieldUnc+self.modelUnc

            #If the Reactor model fission product is not in the Beta models' fission product list
            # save the list of missing branches
            else:
                self.missingCount += thisyield
                self.missingBranch.append(FPZAI)

        # Uncertainty calculation
        #Have to make a 2D Lattice of every single combination of fission products i_j        
        for i in self.FPYlist:
            for j in self.FPYlist:
                #First check that both fission products are in the beta Spectra list
                if i in self.betaSpectraList and j in self.betaSpectraList:
                    adjustmenti = 1
                    adjustmentj = 1
                    #If we're doing an IFP calculation, calculate the adjustments need to be made to 
                    #Both fission products.                    
                    if (ifp_begin < ifp_end):
                        adjustmenti = betaSpectraDB.istplist[i].decay_time_adjust(ifp_begin, ifp_end)
                        adjustmentj = betaSpectraDB.istplist[j].decay_time_adjust(ifp_begin, ifp_end)

                    #Pull the yields, yeild errors, the beta Spectra, and beta Spectral uncertainty
                    #For both fission products (i and j)                    
                    yi = self.FPYlist[i].y*adjustmenti
                    yerri = self.FPYlist[i].yerr*adjustmenti
                    fi = self.betaSpectraList[i]*adjustmenti
                    ferri = self.betaUncertainty[i]*adjustmenti

                    yj = self.FPYlist[j].y*adjustmentj
                    # yerrj = self.FPYlist[j].yerr*adjustmentj
                    fj = self.betaSpectraList[j]*adjustmentj
                    ferrj = self.betaUncertainty[j]*adjustmentj

                    #Carry out the uncertainty calculation                    
                    # if covariance matrix were not loaded, make cov diagonal variance
                    cov_ij = 0
                    if j not in self.FPYlist[i].cov:
                        cov_ij = yerri*yerri if i == j else 0
                    else:
                        cov_ij = self.FPYlist[i].cov[j]

                    sigmay_ij = fi*cov_ij*fj

                    self.uncertainty += sigmay_ij

        # if allowed, add beta model uncertainty to the result
        if modelunc:
            self.uncertainty += self.modelUnc

        self.uncertainty = np.sqrt(self.uncertainty)

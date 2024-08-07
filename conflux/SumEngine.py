# fission product and spectra summation engine

# universal modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
from copy import deepcopy
from tqdm import tqdm

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
    betaUncertainty : dictionary
        A dictionary of the beta Uncertainty's of each isotope
    neutrino : boolean
        A boolean of whether we are looking at the neutrino spectrum or beta spectrum
    xbins : list
        The energy scale of the Summation prediction. Default units are in MeV
    nbins : int
        The number of bins inside your energy scale. More bins means finer calculations
    spectrum : list
        A list to store the total calculated neutrino (beta) spectrum
    uncertainty : list
        A list to store the uncertainties for the neutrino (beta) spectrum

    Methods
    -------

    Clear():
        Clears all dictionaries associated with the SumEngine

    AddModel(fissionModel, W=1.0):
        method to add fission/non-fissile/non-equilibrium isotopes into the engine

    NormalizeFP():
        Normalizes the FPY and the fission product yield Error

    CalcReactorSpectrum(betaSpectraDB, processMissing=False, ifp_begin = 0, ifp_end = 0, modelunc = True):
        Calculates the neutrino Spectrum from the given database with the given energy bins and bounds.
    """

    def __init__(self, neutrino=True, xbins=np.arange(0, 20, 0.1)):
        self.FPYlist = {}
        self.betaSpectraList = {}
        self.betaUncertainty = {}
        self.neutrino = neutrino

        self.xbins = xbins
        self.nbin = len(xbins)
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

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

    # method that accumulates FPYs of fission isotopes in the list of FPY
    def AddFissionIstp(self, isotope, Ei, istpname, IFP=False):
        """
           Add the FPYs of a fission isotope into the the list of FPYs in the model

            Parameters:
                isotope (FissionIstp) : The fission isotope whose products you want to add to the model.
                Ei (float) : The neutron energy that is causing the fissions to occur (0.0, 0.4/0.5, 14)
                fraction (float) : The fractional contribution that this isotope has on the overall model
                d_frac (float) : The uncertainty in the contribution
                IFP (boolean) : determines whether to include the independant fission products in the model

            Returns:
                None
        """

        #Set FPYlist to the Cumulative fission product list of the inputed isotope
        FPYLIST = copy.deepcopy(isotope.CFPY)
        #If we are looking at explosions, or need the independant fission products
        #pull the independant products from the isotope and store it in FPYLIST
        if IFP == True:
            FPYLIST = copy.deepcopy(isotope.IFPY)

        #Check to see if the isotope you want to add has fission products resulting
        #From an interaction with neutrons at the specified energy
        assert Ei in FPYLIST, 'Isotope '+str(isotope.A)+' has only fission data with Ei = ('+str(list(FPYLIST.keys()))+') MeV! (input Ei = ' +str(Ei)+')'

        FPYs = FPYLIST[Ei]
        for FPZAI in FPYs:
            if FPYs[FPZAI].y == 0: continue #If the yield for the specific fission product is 0, skip it (ZAI)

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
            Normalizes the Fission Products in the summation engine.
            Normalization must occur before CalcReactorSpectrum is called
            for the normalization to be applied.

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

    #Calculate the reactor spectrum based off fission yield and beta spectra databases. Options include
    #Processing the missing fission products, calculating the immediate fission products, and including model uncertainties in the uncertainty calculation
    def CalcReactorSpectrum(self, betaSpectraDB, branchErange=[-1, 20], processMissing=False, ifp_begin = 0, ifp_end = 0, modelunc = True, silent = False):
        """
            Calculates the reactor spectrum based off the fission yield database as well as
            the betaSpectra database.
            Parameters:
                betaSpectraDB (BetaEngine): a dictionary of the spectral information for each beta branch
                processMissing (boolean): Determines if missing fission product information is processed
                ifp_begin (float): The start decay time for immediate fission product calculations
                ifp_end (float): The end decay time for immediate fission product calculations
                modelUnc (boolean): Determines if model uncertainties are included in the calculation.
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


        self.betaFPYlist = set(self.FPYlist.keys()).intersection(betaSpectraDB.istplist.keys())
        self.missing_list = set(self.FPYlist.keys()) - set(betaSpectraDB.istplist.keys())

        #Iterate through every fission product in the Reactor Model.
        for FPZAI in tqdm(self.betaFPYlist, desc="Summing beta/neutrino spectra for "+str(len(self.betaFPYlist))+ " fission products", disable=silent):
            # get the yield of the fission products
            thisyield = self.FPYlist[FPZAI].y
            yielderr = self.FPYlist[FPZAI].yerr

            # for IFP calculation, adjust the decay rate of the target isotope
            # by the rate of isotope that are decayed in the time window
            if (ifp_begin < ifp_end):
                adjustment = betaSpectraDB.istplist[FPZAI].decay_time_adjust(ifp_begin, ifp_end)
                thisyield *= adjustment
                yielderr *= adjustment

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

            # #If the Reactor model fission product is not in the Beta models' fission product list
            # # save the list of missing branches
            # else:
            #     self.missingBranch.append(FPZAI)

        # Uncertainty calculation
        #Have to make a 2D Lattice of every single combination of fission products i_j
        for i in tqdm(self.betaFPYlist, desc="Calculating uncertainties of "+str(len(self.betaFPYlist))+ " fission products"):
            for j in self.betaFPYlist:
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

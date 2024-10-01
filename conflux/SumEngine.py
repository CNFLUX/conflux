# fission product and spectra summation engine

# universal modules
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm
import csv

# local modules
from conflux.Basic import Spectrum, Summed
from conflux.BetaEngine import BetaIstp, BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp

class SumEngine(Spectrum, Summed):
    """Class to carry out the summation of reactor antineutrino/beta spectrum by adding the spectra of fission and beta decaying isotopes based on their contributions in the model. The contributions is expected to be provided by the users."""

    FPYlist: dict
    """Dictionary to save fission product yields in the summation model"""
    istplist: dict
    """Dictionary to save all fissile and non-fissile isotope that directly contribute in the summation model"""
    betaSpectraList: dict
    """Dictionary of the spectra of all fissile and non-fissile isotope that directly contribute in the summation model"""
    betaUncertainty: dict
    """Dictionary of the spectrum uncertainty of all fissile and non-fissile isotope that directly contribute in the summation model"""
    countlist: dict
    """Dictionary of contribution (decay accoungting or fractional) of all fissile and non-fissile isotope that directly contribute in the summation model"""
    d_countlist: dict
    """Dictionary of contribution incertainty of all fissile and non-fissile isotope that directly contribute in the summation model"""
    nu_spectrum: bool
    """Whether to calculate neutrino (True) or beta (False) spectrum. default to True"""

    betaDB: BetaEngine
    """Calculated database of beta/neutrino spectra to sum from. The object should only be called after running the :meth:`conflux.BetaEngine.BetaEngine.CalcBetaSpectra` function"""
    
    def __init__(self, betaSpectraDB, nu_spectrum=True):
        """
        Construct the SummationEngine class. The function read an input BetaEngine object that already calculated the beta spectra.
        
        :param betaSpectraDB: A calculated database of beta/neutrino spectra to sum from. The object should only be called after running the :meth:`conflux.BetaEngine.BetaEngine.CalcBetaSpectra` function
        :type betaSpectraDB: :class:`conflux.BetaEngine.BetaEngine`
        :param nu_spectrum: Determine whether to calculate neutrino spectrum, defaults to True
        :type nu_spectrum: bool, optional

        """
        self.FPYlist = {}
        self.istplist = {}
        self.betaSpectraList = {}
        self.modelUncList = {}
        self.yieldUncList = {}
        self.betaUncertainty = {}
        self.countlist = {}
        self.d_countlist = {}

        self.nu_spectrum = nu_spectrum

        self.betaDB = betaSpectraDB
        self.covariance = {}
        
        Spectrum.__init__(self, xbins = self.betaDB.xbins)

    #Self explanatory, clears the various dictionaries associated
    #With the Summation Engine.
    def Clear(self):
        """
        Clear out all associated dictionaries inside the SumEngine for recalculation.
        
        :return: DESCRIPTION
        :rtype: TYPE

        """
        self.FPYlist = {}
        self.istplist = {}
        self.betaUncertainty = {}
        self.modelUncList = {}
        self.yieldUncList = {}
        self.betaSpectraList = {}
        self.countlist = {}
        self.d_countlist = {}
        self.covariance = {}
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

    # method that accumulates FPYs of fission isotopes in the list of FPY
    def AddFissionIstp(self, isotope, istpname, count=1, d_count=0):
        """
        Add the spectrum of a fission isotope into the the list of FPYs in this SumEngine model.
        
        :param isotope: A fission isotope object after running :meth:`conflux.FPYEngine.FissionIstp.CalcBetaSpectra` so the spectrum of the fission isotope can be saved for the summation.
        :type isotope: :class:`conflux.FPYEngine.FissionIstp`
        :param istpname: A unique name assigned to the fission isotope. The name is required to be unique for editting the isotope, instead of redefining the SumEngine model.
        :type istpname: str
        :param count: The counting of total fissions, allowed to be float for fractional calculation, defaults to 1
        :type count: float, optional
        :param d_count: The uncertainty of the counting, defaults to 0
        :type d_count: float, optional

        """
        assert len(isotope.xbins) == self.nbin, "binning of two spectra are different"
        assert istpname not in self.betaSpectraList.keys(), f"{istpname} already exists in the model, please set a different name"

        self.istplist[istpname]=(isotope.id)
        self.betaSpectraList[istpname] = isotope.spectrum
        self.betaUncertainty[istpname] = isotope.uncertainty
        self.yieldUncList[istpname] = isotope.yieldUnc
        self.modelUncList[istpname] = isotope.modelUnc
        self.countlist[istpname] = count
        self.d_countlist[istpname] = d_count
        self.covariance[istpname] = {}
        for key in self.covariance:
            self.covariance[key][istpname] = d_count**2 if key==istpname else 0
            self.covariance[istpname][key] = d_count**2 if key==istpname else 0

        for FPZAI in isotope.FPYlist:
            #Check to see if the fission products in the fission model exists in the Reactor model.

            #If they aren't in the Reactor model, add them to the reactor model, and adjust their
            #Yield and yield error by the weight of the fission mdel.

            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = isotope.FPYlist[FPZAI]
                self.FPYlist[FPZAI].y *= count
                self.FPYlist[FPZAI].yerr *= count


            #If the fission product already exists in the reactor model, take the yield and yield error
            #that this isotope has in this fission model, and add it to the yield and yield error in the
            #Reactor model.
            #Also, Add the covariance of the this fission product from the inputted fission model
            #To the fission product in the Reactor model.
            else:
                self.FPYlist[FPZAI].y += isotope.FPYlist[FPZAI].y*count
                self.FPYlist[FPZAI].yerr += isotope.FPYlist[FPZAI].yerr*count
                self.FPYlist[FPZAI].AddCovariance(isotope.FPYlist[FPZAI])

    def AddBetaIstp(self, betaIstp, istpname, count=1, d_count=0):
        """
        Add the spectrum of a non-fissile beta decaying isotope into the the list of FPYs in this SumEngine model. NOTE: User can only expect the beta/neutrino spectrum of this individual beta decay, no decay chaine following this decay will be added.
        
        :param isotope: A beta isotope object after running :meth:`conflux.BetaEngine.BetaEngine.CalcBetaSpectra` so the spectrum of the fission isotope can be saved for the summation. If the object is a customized beta/neutrino source, one must add and calculate it in self.betaDB before calling it in this function. 
        :type isotope: :class:`conflux.BetaEngine.BetaIstp`
        :param istpname: A unique name assigned to this isotope. The name is required to be unique for editting the isotope, instead of redefining the SumEngine model.
        :type istpname: str
        :param count: The counting of total decay, allowed to be float for fractional calculation, defaults to 1
        :type count: float, optional
        :param d_count: The uncertainty of the counting, defaults to 0
        :type d_count: float, optional

        """
        thisZAI = betaIstp.ZAI
        self.betaDB.istplist[thisZAI] = betaIstp
        self.yieldUncList[istpname] = 0
        self.modelUncList[istpname] = betaIstp.uncertainty
        self.betaSpectraList[istpname] = betaIstp.spectrum
        self.betaUncertainty[istpname] = betaIstp.uncertainty
        self.countlist[istpname] = count
        self.d_countlist[istpname] = d_count
        self.covariance[istpname] = {}
        for key in self.covariance:
            self.covariance[key][istpname] = d_count**2 if key==istpname else 0
            self.covariance[istpname][key] = d_count**2 if key==istpname else 0

    def EditContribution(self, istpname, count, d_count):
        """
        Edit the count and d_count parameters of a fissile isotope or beta isotope that is already added into this SumEngine model.
        
        :param istpname: The UNIQUE name of the isotope assigned when the isotope was added to this SumEngine model.
        :type istpname: str
        :param count: The to-be-updated counting of total decay, allowed to be float for fractional calculation, defaults to 1
        :type count: float, optional
        :param d_count: The uncertainty of the counting, defaults to 0
        :type d_count: float, optional
    

        """
        self.countlist[istpname] = count
        self.d_countlist[istpname] = d_count
        for key in self.covariance:
            self.covariance[key][istpname] = d_count**2 
            self.covariance[istpname][key] = d_count**2 

    def CalcReactorSpectrum(self):
        """Summing all spectra of the added isotopes based on their count and d_count values. Running this function will update the spectrum and uncertainty, if already calculated."""
        #initialize total spectrum & uncertainty
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

        #Initialize model and Yield Uncertainties, a list of missing branches to the total
        #Contribution, the missing yield, and the total yield.
        self.modelUnc = np.zeros(self.nbin)
        self.yieldUnc = np.zeros(self.nbin)

        for istp in self.betaSpectraList.keys():
            counti = self.countlist[istp]
            d_counti = self.d_countlist[istp]
            fi = self.betaSpectraList[istp]
            ferri = self.betaUncertainty[istp]
            merri = self.modelUncList[istp]
            yerri = self.yieldUncList[istp]

            relunc_a = ferri/fi
            relunc_b = d_counti/counti
            fi = self.betaSpectraList[istp]
            ferri = self.betaUncertainty[istp]
            # di = (fi*counti)**2*(reluc_a**2+relunc_b**2)

            self.spectrum += fi*counti
            # self.uncertainty += di
            
            for istp2 in self.betaSpectraList.keys():
                fj = self.betaSpectraList[istp2]
                cov_ij = self.covariance[istp][istp2]

                variance_ij = fi*cov_ij*fj
                
                self.yieldUnc += variance_ij

                if istp==istp2:
                    variance_ij += counti**2*ferri**2
                    self.yieldUnc += counti**2*yerri**2
                    self.modelUnc += counti**2*merri**2
                
                self.uncertainty += variance_ij
        
        self.yieldUnc = np.sqrt(self.yieldUnc)
        self.modelUnc = np.sqrt(self.modelUnc)
        self.uncertainty = np.sqrt(self.uncertainty)
        
    def CustomCovariance(self, DBname, percent = False, rel = False):
        """
           Loads a user defined covariance matrix into the model. The loading of the covariance matrix into the model
           is very similar to the FissionIstp method LoadCovarianceDB

            Parameters:
                DBname (String) : The path to the user defined covariance csv file. has the format "/path/to/file"
                percent (boolean) : Determines whether the covarainces are relative or absolute
                rel (boolean) : Determines whether the product of the yields of each pair of isotopes is a value, or 1
            Returns:
                None
        """

        # determine the whether the covariance are relative or absolute
        rate = 1e4 if percent else 1
        with open(DBname) as inputfile:
            reader = csv.DictReader(inputfile, dialect='excel', delimiter=',')

            for row in reader:
                row_id = int(row[''])
                z = int(row_id/10000)
                i = int(row_id%10)
                a = int((row_id-z*10000-i*1000)/10)
                zai = z*10000+a*10+i
                # keystr = str(key)

                # run through isotopes in saved in the model
                for istpname in self.istplist:
                    if self.istplist[istpname] != zai:
                        continue
                    y1 = self.count(istpname)
                    
                    # find correlated model in the covariance matrix
                    for corristp in self.istplist:    
                        col_ids = self.istplist(corristp)
                        
                        for col_id in col_ids:
                        
                            if str(col_id) not in row:
                                continue 
                            y2 = self.count(corristp)
                            y_prod = (y1*y2)**rel #multiply with y1 and y2 if relative covariance is given
                            self.covariance[istpname][col_id] = y_prod*float(row[str(col_id)])/rate
                            self.covariance[col_id][istpname] = y_prod*float(row[str(istpname)])/rate


    # method to add fission/non-fissile/non-equilibrium isotopes into the engine
    def AddModel(self, fissionModel, W=1.0):
        """
            Adds a reactor model to the summation engine. (obsolete)

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
    def OldCalcReactorSpectrum(self, betaSpectraDB, branchErange=[-1, 20], processMissing=False, ifp_begin = 0, ifp_end = 0, modelunc = True, silent = False):
        """
            Calculates the reactor spectrum based off the fission yield database as well as the betaSpectra database. (obsolete)
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

        # Find common isotopes among beta decaying isotopes and fission products
        self.betaFPYlist = set(self.FPYlist.keys()).intersection(self.betaDB.istplist.keys())
        self.missing_list = set(self.FPYlist.keys()) - set(self.betaDB.istplist.keys())

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

                variance_ij = fi*cov_ij*fj
                self.uncertainty += variance_ij

        self.uncertainty = np.sqrt(self.uncertainty)
        # if allowed, add beta model uncertainty to the result
        if modelunc:
            self.uncertainty += self.modelUnc

# if __name__ == "__main__":
#     # conflux modules
#     from conflux.BetaEngine import BetaEngine
#     from conflux.FPYEngine import FissionModel, FissionIstp

#     xbins = np.arange(0, 20, 0.1)

#     U235 = FissionIstp(92, 235, Ei = 0.5, DB='ENDF', IFPY=False)
#     U235.LoadFissionDB(Ei = 0.5)
#     U235.LoadCorrelation(DB='ENDF')

#     Pu239 = FissionIstp(94, 239, Ei = 0)
#     Pu239.LoadFissionDB()
#     Pu239.LoadCorrelation()

#     model = FissionModel()
#     model.AddContribution(isotope=U235, fraction=1)

#     sum1 = SumEngine(xbins = xbins)
#     sum1.AddModel(model)

#     betaSpectraDB = BetaEngine(xbins=xbins)
#     #betaSpectraDB = BetaEngine(newlist)
#     betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

#     sum1.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False)
#     summed_spect = sum1.spectrum
#     summed_err = sum1.uncertainty
#     summed_model_err = sum1.modelUnc
#     summed_yerr = sum1.yieldUnc

#     fig, ax = plt.subplots()
#     ax.set_xlim([0, 10])
#     ax.set(xlabel='E (MeV)', ylabel='neutrino count')
#     ax.errorbar(sum1.xbins, sum1.spectrum, yerr=sum1.uncertainty, label="test spectrum")
#     ax.legend()
#     plt.show()

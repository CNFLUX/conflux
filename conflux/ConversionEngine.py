# beta spectra conversion engine:
# input: beta spectra measured on nuclear reactors
# output: neutrino spectra of nuclear reactors

# universal modules
import sys
import csv
import numpy as np
from scipy.optimize import curve_fit
from copy import deepcopy
import timeit


# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp

# Class that reads conversion DB and form reference spectra
class BetaData:
    def __init__(self, inputDB, rel_err=True):
        self.x = []
        self.y = []
        self.yerr = []
        self.inputDB = inputDB
        self.LoadConversionDB(inputDB)
        self.spectrum = np.array(self.y)
        self.uncertainty = np.array(self.yerr)

    def LoadConversionDB(self, inputDB, rel_err=True):
        self.x = []
        self.y = []
        self.yerr = []
        with open(inputDB, newline='') as inputCSV:
            inputreader = csv.DictReader(inputCSV, delimiter=',', quotechar='|')
            for row in inputreader:
                E = float(row["E"])
                # convert from keV to MeV
                if E > 100:
                   E /= 1000.

                # convert relative uncertainty in percentage to absolute
                # uncertainty
                y = float(row["Ne"])
                yerr = float(row["dNe"])
                if rel_err:
                    yerr = y*yerr/100.

                self.x.append(E)
                self.y.append(y)
                self.yerr.append(yerr)
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        self.yerr = np.array(self.yerr)

# Class that creates virtual branch based on nuclear data
class VirtualBranch:
    def __init__(self, fisIstp, Ei = 0, Zlist={}, Alist={}, fblist={}, wmlist={}):
        """
        Class that creates virtual branch based on nuclear data
        
        Parameters
        ----------
        fisIstp : FissionIstp
            input fission isotope model to calculate average atom number
        Ei: float
            incedent neutron energy that triggers the fission
        Zlist: dictionary with float keys and float values
            The mapping between the beta energy and input average Z number
        Alist: dictionary with float keys and float values
            The mapping between the beta energy and input average A number
        fblist: dictionary with float keys and float values
            The mapping between the beta energy and customized ratio 1st order
            forbidden transition in the corresponding virtual branches
        wmlist: dictionary with float keys and float values
            The mapping between the beta energy and weak magnatism correction of
            corresponding virtual branches
        """
        
        self.Zavg = 47  # rough avg of Z across all energy
        self.Aavg = 117 # rough avg of A across all energy
        self.Zlist = {}
        self.Alist = {}
        self.fblist = {}
        self.wmlist = {}

        self.fisIstp = fisIstp

        # load FPY of target fission isotope
        self.FPYlist = {}
        if (not self.FPYlist):
            self.LoadFPYList(fisIstp, Ei)

        # load FPY of target fission isotope
        betaEngine = BetaEngine(self.FPYlist)
        
        betaEngine.LoadBetaDB()
        self.istplist = betaEngine.istplist
        
        self._Zlist_cp = deepcopy(Zlist)
        self._Alist_cp = deepcopy(Alist)
        self._fblist_cp = deepcopy(fblist)
        self._wmlist_cp = deepcopy(wmlist)

    # Function to load FPY list
    def LoadFPYList(self, fisIstp, Ei = 0):
        for nuclide in fisIstp.CFPY[Ei]:
            fpNuclide = fisIstp.CFPY[Ei][nuclide]
            if fpNuclide.y == 0: continue
            FPZAI = int(fpNuclide.Z*10000+fpNuclide.A*10+fpNuclide.isomer)

            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = fpNuclide
            else:
                self.FPYlist[FPZAI].y += fpNuclide.y
                self.FPYlist[FPZAI].yerr += fpNuclide.yerr

    # function to precisely calculate average Z, A value of the virtual branch
    def CalcZAavg(self, Elow, Ehigh, missing=False):
        frac_sum = 0
        Afrac_sum = 0
        Zfrac_sum = 0
        for ZAI in self.istplist:
            betaIstp = self.istplist[ZAI]
            if not missing and betaIstp.missing:
                continue
            if (betaIstp.Q >= Elow and betaIstp.Q < Ehigh
                and ZAI in self.FPYlist):
                frac_sum += self.FPYlist[ZAI].y
                Afrac_sum += self.FPYlist[ZAI].y*betaIstp.A
                Zfrac_sum += self.FPYlist[ZAI].y*betaIstp.Z
        Aavg = round(Afrac_sum/frac_sum)
        Zavg = round(Zfrac_sum/frac_sum)
        return Zavg, Aavg

    # define the theoretical beta spectrum shape
    def BetaSpectrum(self, x, E0, contribute, Zavg=None, Aavg=None,
                    nu_spectrum=False, forbiddeness = 0, bAc = 4.7):
        if Zavg is None:
            Zavg = self.Zavg
        if Aavg is None:
            Aavg = self.Aavg
        virtualbata = BetaBranch(Zavg, Aavg, I=0, Q=E0, E0=E0,
                                sigma_E0=0, frac = contribute, sigma_frac = 0,
                                forbiddeness=forbiddeness, bAc=bAc)
        return virtualbata.BetaSpectrum(x, nu_spectrum)*contribute

    # function that fit the reference beta spectrum with virtual brances
    def FitData(self, betadata, slicesize):
        self.contribute = {}
        self.E0 = {}
        self.slicesize = slicesize
        # fill the sub lists as cached slices
        subx = [] # sublist x values
        suby = [] # sublist y values
        subyerr = [] # sublist uncertainty
        xhigh = betadata.x[-1]
        datacache = np.copy(betadata.spectrum) # preserve the data
        for it, x in reversed(list(enumerate(betadata.x))):
            if x < xhigh - slicesize or x == betadata.x[0]:
                subx.append(x)
                suby.append(datacache[it])
                subyerr.append(betadata.uncertainty[it])
                # when the sublist is filled in this slice, do fitting
                if len(subx)>1 and len(suby) == len(subx):
                    if self._Zlist_cp:
                        Zavg = round(np.interp(xhigh-slicesize/2,
                                                list(self._Zlist_cp.keys()),
                                                list(self._Zlist_cp.values())))
                    else:
                        Zavg, _ = self.CalcZAavg(xhigh-slicesize, xhigh)
                    if xhigh not in self.Zlist:
                        self.Zlist[xhigh] = Zavg
                        
                    if self._Alist_cp:
                        Aavg = round(np.interp(xhigh-slicesize/2,
                                                list(self._Alist_cp.keys()),
                                                list(self._Alist_cp.values())))
                    else:
                        _, Aavg = self.CalcZAavg(xhigh-slicesize, xhigh)
                    if xhigh not in self.Alist:
                        self.Alist[xhigh] = Aavg
                        
                    if self._fblist_cp:
                        fbratio = (np.interp(xhigh-slicesize/2,
                                            list(self._fblist_cp.keys()),
                                            list(self._fblist_cp.values())))
                    else:
                        fbratio = 0.0
                    if xhigh not in self.fblist:
                        self.fblist[xhigh] = 0.0
                        
                    if self._wmlist_cp:
                        wm = (np.interp(xhigh-slicesize/2,
                                        list(self._wmlist_cp.keys()),
                                        list(self._wmlist_cp.values())))
                    else:
                        wm = 4.7
                    if xhigh not in self.wmlist:
                        self.wmlist[xhigh] = 4.7

                    # initial guess and boundary setting for parameters
                    tempspec = self.BetaSpectrum(betadata.x, xhigh, 1, Zavg,
                                                Aavg)
                    comparison = (datacache/tempspec)
                    comparison[comparison < 0] = np.inf
                    limit = min(comparison)
                    init_guess = [xhigh, limit/2]
                    fitfunc = (lambda x, e0, c:
                                (1-fbratio)*(self.BetaSpectrum(x, e0, c,
                                                                Zavg=Zavg,
                                                                Aavg=Aavg,
                                                                forbiddeness=0,
                                                                bAc=wm))
                                + fbratio*(self.BetaSpectrum(x, e0, c,
                                                                Zavg=Zavg,
                                                                Aavg=Aavg,
                                                                forbiddeness=1,
                                                                bAc=wm)))
                    
                    popt, pcov = curve_fit(fitfunc, subx, suby,
                                           p0 = init_guess, absolute_sigma=True,
                                           bounds=(0, [xhigh*1.5, limit]))
                    self.contribute[xhigh] = popt[1]
                    self.E0[xhigh] = popt[0]

                    # subtract the best fit spectrum from beta data
                    for i in range(len(betadata.x)):
                        datacache[i] -= self.BetaSpectrum(betadata.x[i],
                                                          popt[0], popt[1])

                    # reset cached slices
                    subx = [] # [x]
                    suby = [] # [datacache[i]]
                    subyerr = [] # [betadata.uncertainty[it]]

                xhigh -= slicesize
            else:
                subx.append(x)
                suby.append(datacache[it])
                subyerr.append(betadata.uncertainty[it])

    # function to calculate summed spectra of virtual branches
    def SumBranches(self, x, thresh = 0, nu_spectrum = True):
        """
        SumBranches(self, x, thresh = 0, nu_spectrum = True)
        
        Function to calculate summed spectra of virtual branches
        
        Parameters
        ----------
        x : ndarray
            The array of x values of the spectrum
        thresh: float
            Overrides the threshold. If threhold > 0, only add branch whose
            energy range is above the threshold
        nu_spectrum: bool
            Overrides the nu_spectrum boolean, detemining whether to sum
            neutrino spectra or beta spectra
        
        Returns
        -------
        result: ndarray
            The summed spectra in the form of ndarray
        """
        result = 0
        for s in self.E0:
            if s > thresh: # if thresh > 0, look at spectra in selected region
                vb = BetaBranch(self.Zlist[s], self.Alist[s],
                                frac=self.contribute[s], I=0, Q = self.E0[s],
                                E0=self.E0[s], sigma_E0=0, sigma_frac=0,
                                forbiddeness=self.fblist[s], bAc=self.wmlist[s])
                result += vb.BetaSpectrum(x, nu_spectrum)*vb.frac
            elif sum(result*x) == 0:
                return result*x
        return result
    
    def Covariance(self, betadata, x, samples = 500, thresh = 0,
                   nu_spectrum = True):
        result = []
        print('Calculating covariance matrix...')
        startTiming = timeit.default_timer()

        for i in range(0, samples):
            toy = deepcopy(betadata)
            for it in range(len(betadata.x)-1, 0, -1):
                toy.y[it] = np.random.normal(toy.y[it], toy.yerr[it])
                
            vbnew = deepcopy(self)
            vbnew.FitData(toy, vbnew.slicesize)
            
            result.append(vbnew.SumBranches(x, thresh, nu_spectrum))
        endTiming = timeit.default_timer()
        runTime = endTiming-startTiming
        print("Finished calculating covairance matrix of "
              + str(samples) + " samples.")
        print("Processing time: "+str(runTime)+" seconds")
        return np.cov(np.transpose(result))

# Class that search for best fit vertual branch and calculate total neutrino
# flux
class ConversionEngine:
    def __init__(self):
        self.betadata = {}
        self.fission_frac = {}
        self.fission_dfrac = {}
        self.fisIstp = {}
        self.vblist = {}

    # load the beta spectrum measurement
    def AddBetaData(self, betadata, fisIstp, name, frac, dfrac=0):
        self.betadata[name] = betadata
        self.fission_frac[name] = frac
        self.fission_dfrac[name] = dfrac
        self.fisIstp[name] = fisIstp

    # Function that lets VB to fit against beta data with user chosen slice size
    def VBfitbeta(self, istp, slicesize = 0.5, Ei = 0, Zlist={}, Alist={}):
        assert istp in self.betadata
        # define the virtual branches to be fit
        vbnew = VirtualBranch(self.fisIstp[istp], Ei, Zlist, Alist)
        vbnew.FitData(self.betadata[istp], slicesize)
        self.vblist[istp] = vbnew
    
    def SummedSpectrum(self, x, nu_spectrum = True, cov_samp = 50):
        """
        SummedSpectrum(self, x, nu_spectrum = True, cov_samp = 50)
        
        Function to sum all best fit spectra calculated in the conversion engine
        
        Parameters
        ----------
        x : ndarray
            The array of x values of the spectrum
        nu_spectrum: bool
            Overrides the nu_spectrum boolean, detemining whether to sum
            neutrino spectra or beta spectra
        cov_samp: int
            The number of MC samples to calculate the covariance matrix. If the
            sample number is below 50, calculation of covariance matrix will be
            skipped.
        
        Returns
        -------
        spectrum: ndarray
            The summed spectrum in the form of ndarray
        uncertainty: ndarray
            The uncertainty array of the summed spectrum
        covariance: ndarray
            The covariance matrix of the of the spectrum on given energy bin
        """
        # Determine whether to calculate the covariance matrix
        calc_cov = cov_samp >= 50
        # Define the function outputs
        nbins_x = len(x)
        spectrum = np.zeros(nbins_x)
        uncertainty = np.zeros(nbins_x)
        covariance = np.zeros((nbins_x, nbins_x))
        
        # Summing all fissile isotopes spectra and uncertainties
        for istp in self.betadata:
            this_frac = self.fission_frac[istp]
            this_vb = self.vblist[istp]
            this_dfrac = self.fission_dfrac[istp]
            
            # Summing spectra with respect to the fraction of fissile isotope
            new_spect = this_vb.SumBranches(x, nu_spectrum)
            spectrum += this_frac*new_spect
        
            # Calculate covariance matrix
            new_cov = (calc_cov
                        *this_vb.Covariance(self.betadata[istp], x=x,
                                            samples=cov_samp,
                                            nu_spectrum=nu_spectrum))
            sqratio1 = new_cov/this_frac**2 if this_frac!=0 else 0
            sqratio2 = np.where(new_spect > 0, this_dfrac**2/new_spect, 0)
            new_cov = sqratio1+np.identity(nbins_x)*sqratio2
            covariance += new_cov
            
        # Get the sqrt diagonal of summed covariance matrix as the uncertainty
        # at each energy bin.
        uncertainty = np.sqrt(covariance.diagonal().copy())
        return spectrum, uncertainty, covariance
    

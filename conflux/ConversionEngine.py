# beta spectra conversion engine:
# input: beta spectra measured on nuclear reactors
# output: neutrino spectra of nuclear reactors

# universal modules
import csv
import numpy as np
from copy import deepcopy
import timeit
import matplotlib.pyplot as plt
from iminuit import Minuit
from scipy.optimize import nnls
import math
from tqdm import tqdm

# local modules
from conflux.Basic import Spectrum
from conflux.config import CONFLUX_DB
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp

# Class that reads conversion DB and form reference spectra
class BetaData(Spectrum):
    """
    A class to save the beta spectrum data to be converted in the computer memory.

    ...

    Attributes
    ----------

    x : list
        A list of the x-axis bins
    y : list
        A list of the y-axis values (bin content)
    yerr : list
        A list of the errors on the y-axis (bin error)
    inputDB : str
        The filename where the spectrum data is saved
    spectrum : :class:`numpy.array`
        An array to store the full beta spectrum
    uncertainty : :class:`numpy.array`
        An array to store the full beta spectrum uncertainty


    Methods
    -------

    LoadConversionDB(self, inputDB, rel_err=True):
        Load the beta spectrum data
    """
    
    def __init__(self, inputDB, rel_err=True):
        """
        Construct the BetaData class.
        
        :param inputDB: The file name of the beta spectrum data to fit against.
        :type inputDB: str
        :param rel_err: Indicate whether the uncertainty of the data is relative (True) or absolute (False), defaults to True
        :type rel_err: bool, optional

        """
        self.x = []
        self.y = []
        self.yerr = []
        self.inputDB = inputDB
        self.LoadConversionDB(inputDB, rel_err)
        self.spectrum = np.array(self.y)
        self.uncertainty = np.array(self.yerr)

    def LoadConversionDB(self, inputDB, rel_err=True):
        """
        Load the beta spectrum data.
        
        :param inputDB: The file name of the beta spectrum data to fit against.
        :type inputDB: str
        :param rel_err: Indicate whether the uncertainty of the data is relative (True) or absolute (False), defaults to True
        :type rel_err: bool, optional

        """
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
        self.spectrum = np.array(self.y)
        self.uncertainty = np.array(self.yerr)

class VirtualBranch(Spectrum):
    """
    A class that creates virtual branch with customized end-point energy and virtual nuclide summarized from nuclides with the similar end-point energy based on nuclear data, and fit the virtual branch against a reference beta spectrum.

    ...

    Attributes
    ----------
    fisIstp : :class:`conflux.FPYEngine.FissionIstp`
        The corresponding fission isotope to fit against. This isotope provides nuclear data to generate a virtual nuclide for the virtual branch.
    FPYlist : dict
        A dictionary of fission product yields of the input fission isotope
    contribute : list
        Best fit contributions of vitual branches
    E0 : list
        Best fit end-point energies of virtual branches
    slicesize : float
        The size the energy slices for virtual branch fitting
    Zlist : dict
        A dictionary of Z numbers of the virtual nuclide calculated for each virtual branch of corresponding energy slice (end-point energy)
    Alist : dict
        A dictionary of A numbers of the virtual nuclide calculated for each virtual branch of corresponding energy slice (end-point energy)
    fblist : dict
        A dictionary of forbidden transition fractions of each virtual branch of corresponding energy slice (end-point energy)
    wmlist : dict
        A dictionary of weak magnetism correction factor of each virtual branch of corresponding energy slice (end-point energy)

    Methods
    -------
    LoadFPYList(self):
        Load the fission products of the input fission isotope
    CalcZAavg(self, Elow, Ehigh, missing=False):
        
    """
    
    def __init__(self, fisIstp, Zlist={}, Alist={}, fblist={}, wmlist={}):
        """
        Construct the VirtualBranch class.
        
        :param fisIstp: The corresponding fission isotope to fit against. This isotope provides nuclear data to generate a virtual nuclide for the virtual branch.
        :type fisIstp: :class:`conflux.FPYEngine.FissionIstp`
        :param Zlist: The mapping between the beta energy and input average Z number, defaults to {}
        :type Zlist: dict, optional
        :param Alist: The mapping between the beta energy and input average A number, defaults to {}
        :type Alist: dict, optional
        :param fblist: The mapping between the beta energy and customized ratio 1st order forbidden transition in the corresponding virtual branches, defaults to {}
        :type fblist: dict, optional
        :param wmlist: The mapping between the beta energy and weak magnatism correction of corresponding virtual branches, defaults to {}
        :type wmlist: dict, optional
        
        """
        self.fisIstp = fisIstp

        self.Zavg = 47  # rough avg of Z across all energy
        self.Aavg = 117 # rough avg of A across all energy
        
        # define empty dictionary of average virtual branch key values to fill 
        # up later
        self.Zlist = {} 
        self.Alist = {}
        self.fblist = {}
        self.wmlist = {}
        
        self._Zlist_cp = deepcopy(Zlist)
        self._Alist_cp = deepcopy(Alist)
        self._fblist_cp = deepcopy(fblist)
        self._wmlist_cp = deepcopy(wmlist)

        # load FPY of target fission isotope
        self.FPYlist = {}
        if (not self.FPYlist):
            self.LoadFPYList()

        # load FPY of target fission isotope
        betaEngine = BetaEngine(self.FPYlist)

        betaEngine.LoadBetaDB()
        self.istplist = betaEngine.istplist

    # Function to load FPY list
    def LoadFPYList(self):
        """Load the fission product yield of the fissile isotope. This is to calculate the average atom numbers of the energy slices."""
        for nuclide in tqdm(self.fisIstp.FPYlist, 
                            desc="Tallying fission products of "+
                                str(self.fisIstp.A)):
            fpNuclide = self.fisIstp.FPYlist[nuclide]
            if fpNuclide.y == 0: continue
            FPZAI = int(fpNuclide.Z*10000 + fpNuclide.A*10 + fpNuclide.isomer)

            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = fpNuclide
            else:
                self.FPYlist[FPZAI].y += fpNuclide.y
                self.FPYlist[FPZAI].yerr += fpNuclide.yerr

    def CalcZAavg(self, Elow, Ehigh, missing=False):
        """
        Calculate the average Z and A number of the slice of the virtual beta spectrum.
        
        :param Elow: lower bondary of the spectrum energy slice
        :type Elow: float
        :param Ehigh: igher bondary of the spectrum energy slice
        :type Ehigh: float
        :param missing: whether to count missing contributions in the beta DB, defaults to False
        :type missing: bool, optional


        """
        frac_sum = 0
        Afrac_sum = 0
        Zfrac_sum = 0
        FB_sum = 0
        for ZAI in self.istplist:
            betaIstp = self.istplist[ZAI]
            if not missing and betaIstp.missing:
                continue
            if (betaIstp.Q >= Elow and betaIstp.Q < Ehigh and
                ZAI in self.FPYlist):
                betaFP = self.FPYlist[ZAI]
                frac_sum += betaFP.y
                Afrac_sum += betaFP.y*betaIstp.A
                Zfrac_sum += betaFP.y*betaIstp.Z
                
                # tallying forbidden branches of a beta decay isotope
                for e0, branch in betaIstp.branches.items():
                    if branch.forbiddenness != 0:
                        FB_sum += branch.frac*betaFP.y
                
        # Calculate the average
        Aavg = round(Afrac_sum/frac_sum)
        Zavg = round(Zfrac_sum/frac_sum)
        FBavg = FB_sum/frac_sum
        return Zavg, Aavg, FBavg

    # define the theoretical beta spectrum shape
    def BetaSpectrum(self, x, E0, contribute, Zavg=None, Aavg=None,
                    nu_spectrum=False, forbiddenness=0, bAc=4.7, norm=True):
        """
        Calculate the theoretical beta/neutrino spectrum shape.
        
        :param x: Energy of the beta/neutrino (x values)
        :type x: float
        :param E0: End-point energy of the spectrum
        :type E0: float
        :param contribute: The normalization or contribution of this spectrum to the summed entire beta spectrum 
        :type contribute: float
        :param Zavg: The Z number of the vitual nuclide for a virtual branch beta spectrum calculation, defaults to None
        :type Zavg: int, optional
        :param Aavg: The Z number of the vitual nuclide for a virtual branch beta spectrum calculation, defaults to None
        :type Aavg: int, optional
        :param nu_spectrum: whether to calculate neutrino (True) or beta (False) spectrum, defaults to False
        :type nu_spectrum: bool, optional
        :param forbiddenness: type of forbidden/allowed transitions , defaults to 0
        :type forbiddenness: int, optional
        :param bAc: Weak magnetism correction factor, defaults to 4.7
        :type bAc: float, optional
        :param norm: Whether to normalize the beta spectrum integral to 1, defaults to True
        :type norm: bool, optional

        """
        if Zavg is None:
            Zavg = self.Zavg
        if Aavg is None:
            Aavg = self.Aavg
        virtualbeta = BetaBranch(Zavg, Aavg, I=0, Q=E0, E0=E0, sigma_E0=0,
                                frac=contribute, sigma_frac=0,
                                forbiddenness=forbiddenness, bAc=bAc)

        binwidths = abs(x[1]-x[0])
        new_spect = virtualbeta.BetaSpectrum(x, nu_spectrum)

        normalize = 1
        full_range = np.arange(0, 20, binwidths)
        full_spect = np.zeros(len(full_range))
        if norm:
            full_spect = virtualbeta.BetaSpectrum(full_range, nu_spectrum)
            normalize = full_spect.sum()

        new_spect /= normalize*binwidths
        return new_spect

    # function that fit the reference beta spectrum with virtual brances
    def FitData(self, betadata, slicesize=0.5):
        """
        Fits the reference beta spectrum with virtual brances data.

        The virtual branches are defined as single beta decay spectra.
        The virtual branches are fitted to slices of the reference data spectra
        in a short range from the highest energy to the lowest, each bestfit
        virtual branch spectrum is subtracted from the reference data to let the
        next virtual branch to fit the updated spectrum in the consequential
        energy slice.
        
        :param betadata: The input BetaData object that provides the reference beta spectrum to fit.
        :type betadata: :class:`conflux.ConversionEngine.BetaData`
        :param slicesize: The size of the energy range, should be greater than 0.2 (MeV) and smaller than 2 (MeV), defaults to 0.5
        :type slicesize: float 

        """
        self.contribute = {}
        self.E0 = {}
        self.slicesize = slicesize
        # fill the sub lists as cached slices
        subx = [] # sublist x values
        suby = [] # sublist y values
        subyerr = [] # sublist uncertainty
        xhigh = 9
        datacache = np.copy(betadata.y) # preserve the data
        fitsequence = reversed(list(enumerate(betadata.x)))
        for it, x in tqdm(fitsequence, 
                          desc="Fitting beta spectrum with reversed sequence"):
            if x <= xhigh - slicesize or x == betadata.x[0]:
                subx.append(x)
                suby.append(datacache[it])
                subyerr.append(betadata.yerr[it])
                # when the sublist is filled in this slice, do fitting
                if len(subx)>1 and len(suby) == len(subx):
                    # setup the virtual isotope properties
                    # find the average Z and A values by interpolating the 
                    # user-given Zlist and Alist
                    if self._Zlist_cp:
                        Zavg = round(np.interp(xhigh-slicesize/2,
                                                list(self._Zlist_cp.keys()),
                                                list(self._Zlist_cp.values())))
                    else:
                        Zavg, _, _ = self.CalcZAavg(xhigh-slicesize, xhigh)
                    if xhigh not in self.Zlist:
                        self.Zlist[xhigh] = Zavg

                    if self._Alist_cp:
                        Aavg = round(np.interp(xhigh-slicesize/2,
                                                list(self._Alist_cp.keys()),
                                                list(self._Alist_cp.values())))
                    else:
                        _, Aavg, _ = self.CalcZAavg(xhigh-slicesize, xhigh)
                    if xhigh not in self.Alist:
                        self.Alist[xhigh] = Aavg

                    if self._fblist_cp:
                        fbratio = (np.interp(xhigh-slicesize/2,
                                            list(self._fblist_cp.keys()),
                                            list(self._fblist_cp.values())))
                    else:
                        _, _, fbratio = self.CalcZAavg(xhigh-slicesize, xhigh)
                        
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

                    subx = np.array(subx)
                    suby = np.array(suby)
                    # suby[suby<0]=0
                    subyerr = np.array(subyerr)

                    # initial guess and boundary setting for parameters
                    # TODO there is a correlation between fitting quality and the range
                    stepsize = 0.02
                    e_upper = xhigh+slicesize
                    e_lower = xhigh-slicesize/8
                    if math.isclose(subx[0], 9, rel_tol = 1e-6):
                        e_upper = xhigh+5

                    betafunc = (lambda x, e:
                                ((1-fbratio)*
                                 (self.BetaSpectrum(x, e, 1, Zavg=Zavg,
                                                    Aavg=Aavg,
                                                    forbiddenness=0,
                                                    bAc=wm))
                                + fbratio*(self.BetaSpectrum(x, e, 1, 
                                                             Zavg=Zavg, 
                                                             Aavg=Aavg,
                                                             forbiddenness=1,
                                                             bAc=wm))))

                    binwidths = betadata.x[1]-betadata.x[0]
                    full_range = np.arange(0, 20, binwidths)

                    best_energy = np.inf
                    best_norm = 0
                    testvalue = np.inf
                    a = 0
                    # find the bestfit virtual branches to the spectrum slices
                    for energy in np.arange(e_lower, e_upper, stepsize):
                        # print(energy)
                        full_spect = betafunc(full_range, energy)
                        norm = full_spect.sum()

                        leastfunc = lambda c: np.sum((suby - 
                                                      c*betafunc(subx, energy)/
                                                      norm/binwidths)**2/
                                                     suby**2)
                        m1 = Minuit(leastfunc, c = 1)

                        m1.migrad()
                        a = m1.values[0]
                        f = m1.fval

                        if (f < testvalue):
                            testvalue = f
                            best_energy = energy
                            if a > 0: best_norm = a


                    self.contribute[xhigh] = best_norm
                    self.E0[xhigh] = best_energy

                    newspect = best_norm*betafunc(betadata.x, best_energy)
                    # least = np.sum((suby - best_norm*betafunc(subx, best_energy))**2)
                    # print(e_lower, e_upper, best_norm, best_energy, least)

                    datacache -= newspect

                    # reset cached slices
                    subx = [] # [x]
                    suby = [] # [datacache[i]]
                    subyerr = [] # [betadata.uncertainty[it]]

                xhigh -= slicesize
            elif x <= xhigh or math.isclose(x, xhigh, rel_tol=1e-6):
                subx.append(x)
                suby.append(datacache[it])
                subyerr.append(betadata.uncertainty[it])

    def FitDataNNLS(self, betadata, slicesize, seeds=2000):
        """
        Fit virtual beta branches to the input beta data using the :meth:`scipy.optimize.nnls` method (maybe obsolete).
        
        :param betadata: The input BetaData object that provides the reference beta spectrum to fit.
        :type betadata: :class:`conflux.ConversionEngine.BetaData`
        :param slicesize: The size of the energy range, should be greater than 0.2 (MeV) and smaller than 2 (MeV), defaults to 0.5
        :type slicesize: float 
        :param seeds: number of MC samples for this calculation, defaults to 2000
        :type seeds: int, optional

        """
        self.contribute = {}
        self.E0 = {}
        self.betadata = betadata
        self.slicesize = slicesize

        # Setting the virtual spectra eneregy ranges
        # To get the x ticks for the beta data in the E range (0, 9]
        xscales = betadata.x[betadata.x<=9.0]
        # To count down from the highest tick to the lowest with the defined
        # stepsize, and let each xhigh be the bondaries of the spectrum slices
        xhigh = (np.arange(xscales[-1], xscales[0], -self.slicesize))

        # To pre-declare the least square value and the bestfit values of each
        # virtual spectrum
        least = np.inf
        bestnorm = np.zeros(len(xhigh))
        beste0 = np.zeros(len(xhigh))

        # setup the virtual isotope properties
        for x in xhigh:
            # Calculate the average Z number of isotopes with end point energy
            # sitting in the specific energy range
            if self._Zlist_cp:
                # if given a list of Z numbers
                Zavg = round(np.interp(x-slicesize/2,
                                        list(self._Zlist_cp.keys()),
                                        list(self._Zlist_cp.values())))
            else:
                Zavg, _ = self.CalcZAavg(x-slicesize, x)
            # fill the Z numbers in the list of virtual isotopes
            if x not in self.Zlist:
                self.Zlist[x] = Zavg

            # calculate the avg A number in the same fashion above
            if self._Alist_cp:
                Aavg = round(np.interp(x-slicesize/2,
                                        list(self._Alist_cp.keys()),
                                        list(self._Alist_cp.values())))
            else:
                _, Aavg = self.CalcZAavg(x-slicesize, x)
            if x not in self.Alist:
                self.Alist[x] = Aavg

            # create the ratio of forbidden transitions in the virtual isotopes
            if self._fblist_cp:
                fbratio = (np.interp(x-slicesize/2,
                                    list(self._fblist_cp.keys()),
                                    list(self._fblist_cp.values())))
            else:
                fbratio = 0.0
            if x not in self.fblist:
                self.fblist[x] = 0.0

            # create the weak magnatism correction fectors for the virtual
            # isotopes
            if self._wmlist_cp:
                wm = (np.interp(x-slicesize/2,
                                list(self._wmlist_cp.keys()),
                                list(self._wmlist_cp.values())))
            else:
                wm = 4.7
            if x not in self.wmlist:
                self.wmlist[x] = 4.7
                
        # Fit the beta spectrum with virtual beta spectra.
        # One virtual spectra for each slice of the beta spectrum.
        # Each beta spectrum's end point energy is randomized to minimize the
        # least square value of the summed virtual spectrum to the beta spectrum
        # data.
        # The virtual branches are fitted with the non-negative least square
        # method to search for the best fit while forbid any below-zero spectrum
        # amplitude.
        bestspectra = []
        for seed in range(seeds):
            # set random end point energy within each energy slice
            randarray = np.random.rand(len(xhigh))
            # limiting the range of parameter randomization
            new_xhigh = xhigh + (randarray*1.)*self.slicesize - 0.3*self.slicesize
            # for the spectrum with highest energy, allow the end point to go upto 12 MeV
            new_xhigh[0] = xhigh[0] + 10*randarray[0]*self.slicesize

            # define the spectra matrix with and fill it virtual spetra
            spectra_matrix = []
            for energy in new_xhigh:
                spectra_matrix.append(self.BetaSpectrum(betadata.x, energy, 1)/betadata.y)

            # find the transpose
            spectra_matrix = np.array(spectra_matrix)
            spectra_matrix = np.transpose(spectra_matrix)
            a, rnorm = nnls(spectra_matrix, betadata.y/betadata.y, maxiter=2000)
            new_spect = (np.dot(spectra_matrix, a))

            # find the minimum within the randomized end point energy groups
            leastfunc = np.sqrt(np.sum((new_spect - 1)**2))
            
            if least>rnorm:
                least = rnorm
                bestnorm = a
                beste0 = new_xhigh

                bestspectra = np.transpose(spectra_matrix)*betadata.y

        self.contribute = dict(zip(xhigh, bestnorm))
        self.E0 = dict(zip(xhigh, beste0))
        self.bestfits =  dict(zip(xhigh, bestspectra))

    # function to calculate summed spectra of virtual branches
    def SumBranches(self, x, thresh = 0, nu_spectrum = False):
        """
        Sum the best fit vitual branch spectra to form the best fit beta/neutrino model spectrum.
        
        :param x: energy of the particle (the x value of spectrum)
        :type x: ndarray
        :param thresh: Overrides the threshold. If threhold > 0, only add branch whose
            energy range is above the threshold, defaults to 0
        :type thresh: float, optional
        :param nu_spectrum: whether to calculate neutrino (True) or beta (False) spectrum, defaults to False
        :type nu_spectrum: bool, optional
        :return: summed spectrum
        :rtype: ndarray

        """
        result = 0

        spect_array = []
        cont_array = []
        # best_array = []
        for s in self.E0:
            if s > thresh: # if thresh > 0, look at spectra in selected region
                vb = BetaBranch(self.Zlist[s], self.Alist[s],
                                frac=1-self.fblist[s], I=0, Q = self.E0[s],
                                E0=self.E0[s], sigma_E0=0, sigma_frac=0,
                                forbiddenness=0, bAc=self.wmlist[s])
                vb_fb = BetaBranch(self.Zlist[s], self.Alist[s],
                                frac=self.fblist[s], I=0, Q = self.E0[s],
                                E0=self.E0[s], sigma_E0=0, sigma_frac=0,
                                forbiddenness=1, bAc=self.wmlist[s])
                newspect = vb.frac * ((1-self.fblist[s])
                                        * vb.BetaSpectrum(x, nu_spectrum)
                                    + self.fblist[s]
                                        * vb_fb.BetaSpectrum(x, nu_spectrum))

                binwidths = x[1]-x[0]
                full_range = np.arange(0, 20, binwidths)
                # this_range = np.arange(x[0], x[-1], 0.01)
                full_spect = vb.frac * ((1-self.fblist[s])
                                    * vb.BetaSpectrum(full_range, nu_spectrum)
                                + self.fblist[s]
                                    * vb_fb.BetaSpectrum(full_range, nu_spectrum))
                # this_spect = vb.frac * ((1-self.fblist[s])
                #                     * vb.BetaSpectrum(this_range, nu_spectrum)
                #                 + self.fblist[s]
                #                     * vb_fb.BetaSpectrum(this_range, nu_spectrum))
                norm = full_spect.sum() if s > binwidths else newspect.sum()
                newspect /= norm*binwidths
                # print('max', max(newspect)/max(self.bestfits[s]))
                # print('sum', sum(newspect)/sum(self.bestfits[s]))

                spect_array.append(newspect)
                cont_array.append(self.contribute[s])
                # best_array.append(self.bestfits[s])
                # fig = plt.figure()
                # plt.title('two best fit compare'+str(s))
                # plt.plot(x, newspect)
                # plt.plot(self.betadata.x, self.bestfits[s])
                # plt.show()

        if len(spect_array) == 0:
            return result*x

        # print(thresh, spect_array, cont_array)
        spect_array = np.transpose(np.array(spect_array))
        cont_array = np.array(cont_array)

        result = np.transpose(spect_array @ cont_array)

        # bestfit = np.transpose(cont_array @ best_array)
        # fig = plt.figure()
        # plt.title('two best fit compare')
        # plt.plot(x, result)
        # plt.plot(self.betadata.x, bestfit)
        # plt.show()


        return result

    def Covariance(self, betadata, x, samples = 50, thresh = 0,
                   nu_spectrum = True):
        """
        Calculate the covariance of the conversion/fitting result by creating MC samples of beta data with varied bin content from the original spectrum within its uncertainty.
        
        :param betadata: The input BetaData object that provides the reference beta spectrum to fit.
        :type betadata: :class:`conflux.ConversionEngine.BetaData`
        :param x: energy of the particle (the x value of spectrum)
        :type x: ndarray
        :param samples: the amount of MC samples, defaults to 50
        :type samples: int, optional
        :param thresh: Overrides the threshold. If threhold > 0, only add branch whose
            energy range is above the threshold, defaults to 0
        :type thresh: float, optional
        :param nu_spectrum: whether to calculate neutrino (True) or beta (False) spectrum, defaults to True
        :type nu_spectrum: bool, optional
        :return: the covariance matrix of the virtual branch fitting or conversion result
        :rtype: ndarray

        """
        result = []
        print(f'Calculating covariance matrix with {samples} MC samples...')
        startTiming = timeit.default_timer()

        # plt.figure()   # TODO comment me
        # vbbuffer = deepcopy(self)
        # vbbuffer.FitData(betadata, vbbuffer.slicesize) # TODO comment me
        # y0 = vbbuffer.SumBranches(x, thresh, nu_spectrum)
        for i in (range(0, samples)):
            toy = deepcopy(betadata)
            for it in range(len(betadata.x)-1, -1, -1):
                toy.y[it] = np.random.normal(toy.y[it], toy.yerr[it])
                if toy.y[it] < 0: toy.y[it] = 0

            vbnew = deepcopy(self)
            vbnew.FitData(toy, vbnew.slicesize)
            y = vbnew.SumBranches(x, thresh, nu_spectrum)
            # plt.plot(x, (y-y0)/y0, color='red') # TODO comment me
            # plt.ylim((-0.2, 0.2))

            result.append(y)
        # plt.show() # TODO comment me
        endTiming = timeit.default_timer()
        runTime = endTiming-startTiming
        print(f"Finished calculating covairance matrix of {samples} samples.")
        print(f"Processing time: {runTime} seconds")
        return np.cov(np.transpose(result))

class ConversionEngine(Spectrum):
    """
    Class that search for best fit vertual branch and sum the converted total neutrino flux.

    ...

    Attributes
    ----------
    betadata : dict
        A dictionary of beta spectra of fission isotopes 
    fission_frac : dict
        A dictionary of the fission isotope fractions
    fission_dfrac : dict
        A dictionary of the fission isotope fraction uncertainties
    fisIstp : dict
        A dictionary of the fission isotope in the conversion model
    vblist: dict
        A dictionary of the virtual branches
    spectrum : :class:`numpy.array`
        An array to store the full beta spectrum
    uncertainty : :class:`numpy.array`
        An array to store the full beta spectrum uncertainty
    covariance : :class:`numpy.array`
        The covariance matrix of the best fit spectrum or converted spectrum

    Methods
    -------
    AddBetaData(self, betadata, fisIstp, name, frac, dfrac=0):
        Add beta spectrum data as fitting reference into the dictionaries
    VBfitbeta(self, istp, slicesize = 0.5, Ei = 0, Zlist={}, Alist={}):
        Let vitual branches to fit against beta data with user chosen slice size.
    SummedSpectrum(self, x, nu_spectrum=True, cov_samp=20):
        Function to sum all best fit spectra calculated in the conversion engine.
    """
    
    def __init__(self):
        """Construct the ConversionEngine class."""
        self.betadata = {}
        self.fission_count = {}
        self.fission_dcount = {}
        self.fisIstp = {}
        self.vblist = {}

    # load the beta spectrum measurement
    def AddBetaData(self, betadata, fisIstp, name, count, d_count=0):
        """
        Add beta spectrum data as fitting reference into the dictionaries.

        :param betadata: The input BetaData object that provides the reference beta spectrum to fit.
        :type betadata: :class:`conflux.ConversionEngine.BetaData`
        :param istp: The corresponding fission isotope to fit against. This isotope provides nuclear data to generate a virtual nuclide for the virtual branch.
        :type istp: :class:`conflux.FPYEngine.FissionIstp`
        :param istp: The unique name of the fission isotope and spectrum that will be fitted.
        :type istp: str
        :param count: The fraction of the fission isotope in the conversion model
        :type count: float
        :param d_count: The fraction uncertainty of the fission isotope in the conversion model
        :type d_count: float

        """
        self.betadata[name] = betadata
        self.fission_count[name] = count
        self.fission_dcount[name] = d_count
        self.fisIstp[name] = fisIstp

    def VBfitbeta(self, istp, slicesize = 0.5, Zlist={}, Alist={}):
        """
        Let vitual branches to fit against beta data with user chosen slice size.
            
        :param istp: The corresponding fission isotope to fit against. This isotope provides nuclear data to generate a virtual nuclide for the virtual branch.
        :type istp: :class:`conflux.FPYEngine.FissionIstp`
        :param slicesize: The size of the energy range, should be greater than 0.2 (MeV) and smaller than 2 (MeV), defaults to 0.5
        :type slicesize: float 
        :param Zlist: The mapping between the beta energy and input average Z number, defaults to {}
        :type Zlist: dict, optional
        :param Alist: The mapping between the beta energy and input average A number, defaults to {}
        :type Alist: dict, optional

        """
        assert istp in self.betadata
        # define the virtual branches to be fit
        vbnew = VirtualBranch(self.fisIstp[istp], Zlist, Alist)
        vbnew.FitData(self.betadata[istp], slicesize)
        self.vblist[istp] = vbnew

    def SummedSpectrum(self, x, nu_spectrum=True, cov_samp=50):
        """
        Sum all best fit  or converted spectra calculated in the conversion engine based on the fractions and uncertainties.
        
        :param x: Particle energy (x values of the specrrum)
        :type x: ndarray
        :param nu_spectrum: whether to calculate beta (False) or netutrino (True) spectrum, defaults to True
        :type nu_spectrum: bool, optional
        :param cov_samp: the amount of MC samples for covariance matrix calculation, defaults to 50
        :type cov_samp: int, optional
        
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
        calc_cov = cov_samp >= 20
        # Define the function outputs
        nbins_x = len(x)
        self.spectrum = np.zeros(nbins_x)
        self.uncertainty = np.zeros(nbins_x)
        self.covariance = np.zeros((nbins_x, nbins_x))

        # Summing all fissile isotopes spectra and uncertainties by looping 
        # through all virtual branches
        for istp in self.betadata:
            this_vb = self.vblist[istp] 
            this_count = self.fission_count[istp]
            this_dcount = self.fission_dcount[istp]

            # Summing spectra with respect to the fraction of fissile isotope
            new_spect = this_vb.SumBranches(x, nu_spectrum=nu_spectrum)
            self.spectrum += this_count*new_spect

            # Calculate covariance matrix
            new_cov = np.zeros((nbins_x, nbins_x))
            if calc_cov:
                new_cov = (this_vb.Covariance(self.betadata[istp], x=x,
                                                samples=cov_samp,
                                                nu_spectrum=nu_spectrum))
            # Propagating the uncertainty of fraction and spectra
            sqratio1 = new_cov/this_count**2 if this_count!=0 else 0
            sqratio2 = np.where(new_spect>0, this_dcount**2/new_spect**2, 0)
            new_cov = sqratio1 + np.identity(nbins_x) * sqratio2
            self.covariance += new_cov

        # Get the sqrt diagonal of summed covariance matrix as the uncertainty
        # at each energy bin.
        self.uncertainty = np.sqrt(self.covariance.diagonal().copy())
        return self.spectrum, self.uncertainty, self.covariance

def HuberZavg(x, c0, c1, c2):
    return c0+c1*x+c2*x**2

def rebin(data, bins):
    # Set the number of points in the averaged array
    num_points = 10

    # Calculate the size of each chunk
    chunk_size = len(data) // num_points

    # Reshape the data into chunks and calculate the average for each chunk
    averaged_data = np.mean(data[:num_points * 
                                 chunk_size].reshape((num_points, -1)), 
                            axis=1)

    return averaged_data

# # testing script
# if __name__ == "__main__":
#     from conflux.FPYEngine import FissionModel, FissionIstp
#     from conflux.SumEngine import SumEngine
#     # from conflux.ConversionEngine import ConversionEngine, BetaData
#     import matplotlib.pyplot as plt


#     # Begin the calculation by sourcing the default beta data
#     beta235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")
#     # beta2351 = BetaData("./data/conversionDB/Synthetic_235_beta.csv")
#     beta235s = BetaData(CONFLUX_DB+"/example_models/U235_synth_data_1.5_9.6.csv")
#     beta239 = BetaData(CONFLUX_DB+"/conversionDB/Pu_239_e_2014.csv")
#     beta241 = BetaData(CONFLUX_DB+"/conversionDB/Pu_241_e_2014.csv")

#     # Define isotopic fission yield DB to calculate average atom numbers of
#     # virtual branches
#     U235 = FissionIstp(92, 235, Ei=0)
#     Pu239 = FissionIstp(94, 239, Ei=0)
#     Pu241 = FissionIstp(94, 241, Ei=0)

#     # Loading default fission product DB
#     U235.LoadFissionDB(DB='JEFF')
#     Pu239.LoadFissionDB()
#     Pu241.LoadFissionDB()

#     # Define the size of energy slice
#     branch_slice = 0.25
#     # Declare the conversion engine by adding beta data with corresponding FPY
#     # database
#     convertmodel = ConversionEngine()
#     # Add beta spectra and fission products to the conversion engine
#     convertmodel.AddBetaData(beta235, U235, "U235", 1.0)
#     convertmodel.VBfitbeta("U235", branch_slice)

#     xval = np.arange(0, 10, 0.01)
#     #Something wrong with the plotting on this one
#     for i in range(0, 20):
#         print(i*0.5, 
#               sum(convertmodel.vblist["U235"].SumBranches(xval, 
#                                                           thresh =i*0.5, 
#                                                           nu_spectrum = False)))
#         if not sum(convertmodel.vblist["U235"].SumBranches(xval, 
#                                                            thresh =i*0.5, 
#                                                            nu_spectrum = False))>=0:
#             print("sth wrong i die")
#             continue
#         plt.errorbar(xval, 
#                      convertmodel.vblist["U235"].SumBranches(xval, 
#                                                             thresh =i*0.5, 
#                                                             nu_spectrum = False), 
#                      fmt='--', label='test')
#     #Plot out the raw beta data
#     plt.errorbar(convertmodel.betadata["U235"].x, 
#                  convertmodel.betadata["U235"].y, 
#                  convertmodel.betadata["U235"].yerr, label='beta data')
#     #Plot out the calculated beta spectrum
#     plt.errorbar(xval, 
#                  convertmodel.vblist["U235"].SumBranches(xval, 
#                                                          nu_spectrum = False), 
#                  label='beta')
#     #Plot out the calculated neutrino spectrum
#     plt.errorbar(xval, 
#                  convertmodel.vblist["U235"].SumBranches(xval, 
#                                                          nu_spectrum = True), 
#                  label='neutrino')
#     plt.legend()
#     plt.show()

#     Zlist = dict(zip(xval, HuberZavg(xval, 49, -0.4, -0.084)))
#     #convertmodel.VBfitbeta("U235", branch_slice)

#     final_spect, final_unc, final_cov = convertmodel.SummedSpectrum(xval, 
#                                                                     nu_spectrum=False, 
#                                                                     cov_samp=20)
#     final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval, 
#                                                                        nu_spectrum=True, 
#                                                                        cov_samp=20)

#     newxval = np.arange(2.125, 8.375, 0.25)
#     # newyval = Rebin(xval, final_spect, newxval)
#     newyval1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(newxval, 
#                                                                    nu_spectrum=True, 
#                                                                    cov_samp=5)
#     # for i in convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True):
#     #     print(i)

#     # testxval = np.linspace(0,200,201)
#     # testyval = 1*testxval+2
#     # testxout = np.linspace(0,200,21)
#     # testyout = Rebin(testxval, testyval, testxout)
#     # for a in (newyval1):
#     #     print(a)

#     # fig = plt.figure()
#     # # plt.yscale('log')
#     # plt.errorbar(convertmodel.betadata["U235"].x,
#     #     convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
#     #     label='beta data')
#     # plt.plot(newxval, newyval1, label='neutrino rebined')
#     # plt.plot(xval, final_spect1, label='best fit neutrino')
#     # plt.legend()
#     # plt.show()
#     # fig.savefig("bestfit_spectra.png")

#     betaspect = np.interp(xval, 
#                           convertmodel.betadata["U235"].x, 
#                           convertmodel.betadata["U235"].y)
#     diff = (final_spect-betaspect)/betaspect
#     fig = plt.figure()
#     plt.title('Beta fit beta residual')
#     plt.ylim([-0.1, 0.1])
#     plt.xlim([2,9])
#     plt.plot(xval, diff)
#     # plt.plot(xval, betaspect)
#     # plt.plot(convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y, 'o')

#     plt.xlabel('Beta E (MeV)')
#     plt.ylabel('Residual')
#     plt.show()
#     fig.savefig("bestfit_beta_compare.png")

#     betaspect = np.interp(xval, 
#                           convertmodel.betadata["U235"].x, 
#                           convertmodel.betadata["U235"].y)
#     diff = (final_spect-betaspect)/betaspect
#     print('leastsquare', sum(diff[xval<=9]**2))
#     fig = plt.figure()
#     plt.title('Beta and neutrino spectra')
#     plt.yscale('log')
#     plt.xlim([0, 10])
#     plt.plot(xval, final_spect, label='bestfit beta')
#     plt.plot(xval, betaspect, label='beta data')
#     plt.plot(convertmodel.betadata["U235"].x, 
#              convertmodel.betadata["U235"].y, 
#              'o')

#     plt.xlabel('Beta E (MeV)')
#     plt.ylabel('Residual')
#     plt.legend()
#     plt.show()
#     fig.savefig("bestfit_beta_compare2.png")
#     neu235s = BetaData(CONFLUX_DB+"/example_models/U235_synth_compare.csv")

#     new235spect = neu235s.y
#     diff = (final_spect1-new235spect)/new235spect
#     final_spect1_avg = rebin(final_spect1, 40)
#     new235spect_avg = rebin(new235spect, 40)
#     xval_avg = rebin(xval, 40)
#     diff_avg = (final_spect1_avg-new235spect_avg)/new235spect_avg
#     fig = plt.figure()
#     plt.title("converted neutrino compared to synthetic neutrino spectrum")
#     plt.ylim([-0.2, 0.2])
#     plt.xlim([2,9])
#     plt.plot(xval, diff)
#     plt.plot(xval_avg, diff_avg)
#     plt.xlabel('E (MeV)')
#     plt.ylabel('Residual')
#     plt.show()
#     fig.savefig("bestfit_neu_compare.png")

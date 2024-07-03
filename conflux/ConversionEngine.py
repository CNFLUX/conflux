# beta spectra conversion engine:
# input: beta spectra measured on nuclear reactors
# output: neutrino spectra of nuclear reactors

# universal modules
import sys
import csv
import numpy as np
from scipy.optimize import curve_fit, nnls
from scipy.stats import chisquare
from scipy import interpolate
from copy import deepcopy
import timeit
import matplotlib.pyplot as plt
from iminuit.cost import LeastSquares
from iminuit import Minuit
import math
from tqdm import tqdm

# local modules
from conflux.config import CONFLUX_DB
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
        self.spectrum = np.array(self.y)
        self.uncertainty = np.array(self.yerr)

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
        for nuclide in tqdm(self.fisIstp.CFPY[Ei], desc="Tallying fission products of "+str(self.fisIstp.A)):
            fpNuclide = fisIstp.CFPY[Ei][nuclide]
            if fpNuclide.y == 0: continue
            FPZAI = int(fpNuclide.Z*10000 + fpNuclide.A*10 + fpNuclide.isomer)

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
                    nu_spectrum=False, forbiddenness=0, bAc=4.7, norm=True):
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
            #normalize = new_spect.sum()*full_spect.sum()/this_spect.sum() if E0 > binwidths else new_spect.sum()
            # print("test", this_spect.sum()/full_spect.sum())

        new_spect /= normalize*binwidths
        # print("again", sum(new_spect))
        return new_spect

    # function that fit the reference beta spectrum with virtual brances
    def FitData(self, betadata, slicesize):
        """
        Fits the reference beta spectrum with virtual brances data.

        The virtual branches are defined as single beta decay spectra.
        The virtual branches are fitted to slices of the reference data spectra
        in a short range from the highest energy to the lowest, each bestfit
        virtual branch spectrum is subtracted from the reference data to let the
        next virtual branch to fit the updated spectrum in the consequential
        energy slice.

        Parameters
        ----------
        betadata : array
        slicesize : float
            The size of the energy range, should be greater than 0.2 (MeV) and
            smaller than 2 (MeV)
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
        for it, x in tqdm(fitsequence, desc="Fitting beta spectrum with reversed sequence"):
            if x <= xhigh - slicesize or x == betadata.x[0]:
                subx.append(x)
                suby.append(datacache[it])
                subyerr.append(betadata.yerr[it])
                # when the sublist is filled in this slice, do fitting
                if len(subx)>1 and len(suby) == len(subx):
                    # setup the virtual isotope properties
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
                                ((1-fbratio)*(self.BetaSpectrum(x, e, 1,
                                                                Zavg=Zavg,
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

                        leastfunc = lambda c: np.sum((suby - c*betafunc(subx, energy)/norm/binwidths)**2/suby**2)
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
                    least = np.sum((suby - best_norm*betafunc(subx, best_energy))**2)
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

    # TODO: revisit the NNLS fitter by changing the selection strat of the energy slice
    def FitDataNNLS(self, betadata, slicesize, seeds=2000):
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

            # TODO: There is one problem, the least square found by NNLS does
            # not produce least value in a separated calculation
            if least>rnorm:
                least = rnorm
                bestnorm = a
                beste0 = new_xhigh

                bestspectra = np.transpose(spectra_matrix)*betadata.y

        self.contribute = dict(zip(xhigh, bestnorm))
        self.E0 = dict(zip(xhigh, beste0))
        self.bestfits =  dict(zip(xhigh, bestspectra))
        # subtract the best fit spectrum from final_spect
        fig = plt.figure()
        # plt.yscale('log')
        # plt.ylim([-1.5*max(abs(suby)), 1.5*max(abs(suby))])
        # plt.errorbar(betadata.x, betadata.spectrum, betadata.uncertainty)
        plt.plot(betadata.x, new_spect - 1)

        plt.pause(1)
        plt.show()
        print(least, self.E0, self.contribute)

    # function to calculate summed spectra of virtual branches
    def SumBranches(self, x, thresh = 0, nu_spectrum = False):
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

        spect_array = []
        cont_array = []
        best_array = []
        for s in self.E0:
            if s > thresh: # if thresh > 0, look at spectra in selected region
                vb = BetaBranch(self.Zlist[s], self.Alist[s],
                                frac=1, I=0, Q = self.E0[s],
                                E0=self.E0[s], sigma_E0=0, sigma_frac=0,
                                forbiddenness=0, bAc=self.wmlist[s])
                vb_fb = BetaBranch(self.Zlist[s], self.Alist[s],
                                frac=1, I=0, Q = self.E0[s],
                                E0=self.E0[s], sigma_E0=0, sigma_frac=0,
                                forbiddenness=1, bAc=self.wmlist[s])
                newspect = vb.frac * ((1-self.fblist[s])
                                        * vb.BetaSpectrum(x, nu_spectrum)
                                    + self.fblist[s]
                                        * vb_fb.BetaSpectrum(x, nu_spectrum))

                binwidths = x[1]-x[0]
                full_range = np.arange(0, 20, binwidths)
                this_range = np.arange(x[0], x[-1], 0.01)
                full_spect = vb.frac * ((1-self.fblist[s])
                                        * vb.BetaSpectrum(full_range, nu_spectrum)
                                    + self.fblist[s]
                                        * vb_fb.BetaSpectrum(full_range, nu_spectrum))
                this_spect = vb.frac * ((1-self.fblist[s])
                                        * vb.BetaSpectrum(this_range, nu_spectrum)
                                    + self.fblist[s]
                                        * vb_fb.BetaSpectrum(this_range, nu_spectrum))
                norm = full_spect.sum() if s > binwidths else newspect.sum()
                # print('e0', s, full_spect.sum()/this_spect.sum(), newspect.sum(), sum(newspect))
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

            elif sum(result*x) == 0:
                return result*x

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

    def Covariance(self, betadata, x, samples = 500, thresh = 0,
                   nu_spectrum = True):
        result = []
        print('Calculating covariance matrix...')
        startTiming = timeit.default_timer()

        for i in range(0, samples):
            toy = deepcopy(betadata)
            for it in range(len(betadata.x)-1, -1, -1):
                toy.y[it] = np.random.normal(toy.y[it], toy.yerr[it])
                if toy.y[it] < 0: toy.y[it] = 0

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
        # vbnew.FitData(self.betadata[istp], slicesize)
        self.vblist[istp] = vbnew

    def SummedSpectrum(self, x, nu_spectrum=True, cov_samp=50):
        """
        SummedSpectrum(self, x, nu_spectrum=True, cov_samp=50)

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
            new_spect = this_vb.SumBranches(x, nu_spectrum=nu_spectrum)
            spectrum += this_frac*new_spect

            # Calculate covariance matrix
            new_cov = np.zeros((nbins_x, nbins_x))
            if calc_cov:
                new_cov = (calc_cov
                            *this_vb.Covariance(self.betadata[istp], x=x,
                                                samples=cov_samp,
                                                nu_spectrum=nu_spectrum))
            sqratio1 = new_cov/this_frac**2 if this_frac!=0 else 0
            sqratio2 = np.where(new_spect>0, this_dfrac**2/new_spect, 0)
            new_cov = sqratio1 + np.identity(nbins_x) * sqratio2
            covariance += new_cov

        # Get the sqrt diagonal of summed covariance matrix as the uncertainty
        # at each energy bin.
        uncertainty = np.sqrt(covariance.diagonal().copy())
        return spectrum, uncertainty, covariance

def HuberZavg(x, c0, c1, c2):
    return c0+c1*x+c2*x**2

def rebin(data, bins):
    # Set the number of points in the averaged array
    num_points = 10

    # Calculate the size of each chunk
    chunk_size = len(data) // num_points

    # Reshape the data into chunks and calculate the average for each chunk
    averaged_data = np.mean(data[:num_points * chunk_size].reshape((num_points, -1)), axis=1)

    return averaged_data

# testing script
if __name__ == "__main__":
    from conflux.BetaEngine import BetaEngine, BetaBranch
    from conflux.FPYEngine import FissionModel, FissionIstp
    from conflux.SumEngine import SumEngine
    # from conflux.ConversionEngine import ConversionEngine, BetaData
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    # Begin the calculation by sourcing the default beta data
    beta235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")
    # beta2351 = BetaData("./data/conversionDB/Synthetic_235_beta.csv")
    beta235s = BetaData(CONFLUX_DB+"/conflux/conflux/U235_synth_data_1.5_9.6.csv")
    beta239 = BetaData(CONFLUX_DB+"/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData(CONFLUX_DB+"/conversionDB/Pu_241_e_2014.csv")

    # Define isotopic fission yield DB to calculate average atom numbers of
    # virtual branches
    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)

    # Loading default fission product DB
    U235.LoadFissionDB(defaultDB='JEFF')
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()

    # Define the size of energy slice
    branch_slice = 0.25
    # Declare the conversion engine by adding beta data with corresponding FPY
    # database
    convertmodel = ConversionEngine()
    # Add beta spectra and fission products to the conversion engine
    convertmodel.AddBetaData(beta235s, U235, "U235", 1.0)
    print(beta235s)
    # convertmodel.AddBetaData(beta239, Pu239, "Pu239", 1.0)
    # convertmodel.AddBetaData(beta241, Pu241, "Pu241", 1.0)
    # Do virtual branch fitting with the defined virtual branch energy range
    xval = np.arange(0,10, 0.01)
    Zlist = dict(zip(xval, HuberZavg(xval, 49, -0.4, -0.084)))
    #convertmodel.VBfitbeta("U235", branch_slice)
    convertmodel.VBfitbeta("U235", branch_slice)

    final_spect, final_unc, final_cov = convertmodel.SummedSpectrum(xval, nu_spectrum=False, cov_samp=5)
    final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval, nu_spectrum=True, cov_samp=5)

    newxval = np.arange(2.125, 8.375, 0.25)
    # newyval = Rebin(xval, final_spect, newxval)
    newyval1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(newxval, nu_spectrum=True, cov_samp=5)
    # for i in convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True):
    #     print(i)

    # testxval = np.linspace(0,200,201)
    # testyval = 1*testxval+2
    # testxout = np.linspace(0,200,21)
    # testyout = Rebin(testxval, testyval, testxout)
    # for a in (newyval1):
    #     print(a)

    # fig = plt.figure()
    # # plt.yscale('log')
    # plt.errorbar(convertmodel.betadata["U235"].x,
    #     convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
    #     label='beta data')
    # plt.plot(newxval, newyval1, label='neutrino rebined')
    # plt.plot(xval, final_spect1, label='best fit neutrino')
    # plt.legend()
    # plt.show()
    # fig.savefig("bestfit_spectra.png")

    betaspect = np.interp(xval, convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y)
    diff = (final_spect-betaspect)/betaspect
    fig = plt.figure()
    plt.title('Beta fit beta residual')
    plt.ylim([-0.1, 0.1])
    plt.xlim([2,9])
    plt.plot(xval, diff)
    # plt.plot(xval, betaspect)
    # plt.plot(convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y, 'o')

    plt.xlabel('Beta E (MeV)')
    plt.ylabel('Residual')
    plt.show()
    fig.savefig("bestfit_beta_compare.png")

    betaspect = np.interp(xval, convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y)
    diff = (final_spect-betaspect)/betaspect
    print('leastsquare', sum(diff[xval<=9]**2))
    fig = plt.figure()
    plt.title('Beta and neutrino spectra')
    plt.yscale('log')
    plt.xlim([0, 10])
    plt.plot(xval, final_spect, label='bestfit beta')
    plt.plot(xval, betaspect, label='beta data')
    plt.plot(convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y, 'o')

    plt.xlabel('Beta E (MeV)')
    plt.ylabel('Residual')
    plt.legend()
    plt.show()
    fig.savefig("bestfit_beta_compare2.png")
    neu235s = BetaData(os.environ["HOME"]+"/conflux/conflux/examples/U235_synth_compare.csv")

    new235spect = neu235s.y
    diff = (final_spect1-new235spect)/new235spect
    final_spect1_avg = rebin(final_spect1, 40)
    new235spect_avg = rebin(new235spect, 40)
    xval_avg = rebin(xval, 40)
    diff_avg = (final_spect1_avg-new235spect_avg)/new235spect_avg
    fig = plt.figure()
    plt.title("converted neutrino compared to synthetic neutrino spectrum")
    plt.ylim([-0.2, 0.2])
    plt.xlim([2,9])
    plt.plot(xval, diff)
    plt.plot(xval_avg, diff_avg)
    plt.xlabel('E (MeV)')
    plt.ylabel('Residual')
    plt.show()
    fig.savefig("bestfit_neu_compare.png")

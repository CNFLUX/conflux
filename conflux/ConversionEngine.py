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

                # convert relative uncertainty in percentage to absolute uncertainty
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
    def __init__(self, fisIstp, Ei = 0):
        self.Zavg = 47  # rough avg of Z across all energy
        self.Aavg = 117 # rough avg of A across all energy
        self.fisIstp = fisIstp

        # load FPY of target fission isotope
        self.FPYlist = {}
        if (not self.FPYlist):
            self.LoadFPYList(fisIstp, Ei)

        # load FPY of target fission isotope
        betaEngine = BetaEngine(self.FPYlist)
        
        betaEngine.LoadBetaDB()
        self.istplist = betaEngine.istplist

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
            if betaIstp.Q >= Elow and betaIstp.Q < Ehigh:
                frac_sum += self.FPYlist[ZAI].y
                Afrac_sum += self.FPYlist[ZAI].y*betaIstp.A
                Zfrac_sum += self.FPYlist[ZAI].y*betaIstp.Z
        self.Aavg = Afrac_sum/frac_sum
        self.Zavg = Zfrac_sum/frac_sum

    # define the theoretical beta spectrum shape
    def BetaSpectrum(self, x, E0, contribute, nu_spectrum=False, forbiddeness = 0, WM = 0.0047):
        virtualbata = BetaBranch(self.Zavg, self.Aavg, I=0, Q=E0, E0=E0, sigma_E0=0, frac = contribute, sigma_frac = 0, forbiddeness=forbiddeness, WM=WM)
        return virtualbata.BetaSpectrum(x, nu_spectrum)*contribute

    # function that fit the reference beta spectrum with virtual brances
    def FitData(self, betadata, slicesize):
        self.contribute = {}
        self.E0 = {}
        self.Zlist = {}
        self.Alist = {}
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
                if len(subx)>1 and len(suby) == len(subx) :
                    self.CalcZAavg(xhigh-slicesize, xhigh)
                    self.Zlist[xhigh] = self.Zavg
                    self.Alist[xhigh] = self.Aavg

                    # initial guess and boundary setting for parameters
                    tempspec = self.BetaSpectrum(betadata.x, xhigh, 1)
                    comparison = (datacache/tempspec)
                    comparison[comparison < 0] = np.inf
                    limit = min(comparison)
                    init_guess = [xhigh, limit/2]
                    
                    popt, pcov = curve_fit(self.BetaSpectrum, subx, suby, p0 = init_guess, absolute_sigma=True, bounds=(0, [xhigh*1.5, limit]))
                    self.contribute[xhigh] = popt[1]
                    self.E0[xhigh] = popt[0]

                    # subtract the best fit spectrum from beta data
                    for i in range(len(betadata.x)):
                        datacache[i] -= self.BetaSpectrum(betadata.x[i], popt[0], popt[1])

                    # reset cached slices
                    subx = []
                    suby = []
                    subyerr = []

                xhigh -= slicesize
            else:
                subx.append(x)
                suby.append(datacache[it])
                subyerr.append(betadata.uncertainty[it])

    # function to calculate summed spectra of virtual branches
    def SumBranches(self, x, thresh = 0, nu_spectrum = True):
        result = 0
        for s in self.E0:
            if s > thresh: # if thresh > 0, look at spectra in selected region
                vb = BetaBranch(self.Zlist[s], self.Alist[s], frac=self.contribute[s], I=0, Q = self.E0[s], E0=self.E0[s], sigma_E0=0, sigma_frac=0, forbiddeness=0, WM=0.0047)
                result += vb.BetaSpectrum(x, nu_spectrum)*vb.frac
            elif sum(result*x) == 0:
                return result*x
        return result
    
    def Covariance(self, betadata, x, samples = 1000, thresh = 0, nu_spectrum = True):
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
        print("Finished calculating covairance matrix of "+ str(samples) + " samples.")
        print("Processing time: "+str(runTime)+" seconds")
        return np.cov(np.transpose(result))

# class that search for best fit vertual branch and calculate total neutrino flux
class ConversionEngine:
    def __init__(self):
        self.betadata = {}
        self.fissionfrac = {}
        self.fisIstp = {}
        self.vblist = {}

    # load the beta spectrum measurement
    def AddBetaData(self, betadata, fisIstp, name, frac):
        self.betadata[name] = betadata
        self.fissionfrac[name] = frac
        self.fisIstp[name] = fisIstp

    # Function that lets VB to fit against beta data with user chosen slice size
    def VBfit(self, slicesize = 0.5):
        self.slicesize = slicesize
        for istp in self.betadata:
            # define the virtual branches to be fit
            vbnew = VirtualBranch(self.fisIstp[istp])
            vbnew.FitData(self.betadata[istp], slicesize)
            self.vblist[istp] = vbnew
    
        
    

# beta spectra conversion engine:
# input: beta spectra measured on nuclear reactors
# output: neutrino spectra of nuclear reactors

# universal modules
import sys
import csv
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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
        self.LoadDB(inputDB)

    def LoadDB(self, inputDB, rel_err=True):
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

        # load FPY of target fission isotope
        self.FPYlist = {}
        self.LoadFPYList(fisIstp, Ei)

        # load FPY of target fission isotope
        betaSpectra = BetaEngine(self.FPYlist)
        betaSpectra.LoadBetaDB()
        self.betaIstpList = betaSpectra.istplist

    # Function to load FPY list
    def LoadFPYList(self, fisIstp, Ei = 0):
        for nuclide in fisIstp.CFPY[Ei]:
            if nuclide.y == 0: continue
            FPZAI = int(nuclide.Z*10000+nuclide.A*10+nuclide.isomer)

            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = nuclide
            else:
                self.FPYlist[FPZAI].y += nuclide.y
                self.FPYlist[FPZAI].yerr += nuclide.yerr

    # function to precisely calculate average Z, A value of the virtual branch
    def CalcZAavg(self, Elow, Ehigh):
        frac_sum = 0
        Afrac_sum = 0
        Zfrac_sum = 0
        for ZAI in self.betaIstpList:
            for branch in self.betaIstpList[ZAI]:
                if branch.E0 >= Elow and branch.E0 <Ehigh:
                    frac_sum += branch.frac
                    Afrac_sum += branch.frac*branch.A
                    Zfrac_sum += branch.frac*branch.Z
        self.Aavg = Afrac_sum/frac_sum
        self.Zavg = Zfrac_sum/frac_sum

    # define the theoretical beta spectrum shape
    def BetaSpectrum(self, x, E0, contribute, forbiddeness = 0, WM = 0.0047):
        virtualbata = BetaBranch(self.Zavg, self.Aavg, frac=contribute, I=0, E0=E0, sigma_E0=0, forbiddeness=forbiddeness, WM=WM)
        return virtualbata.BetaSpectrum(x)

    # define the theoretical neutrino spectrum shape
    def NueSpectrum(self, x, E0, contribute, forbiddeness = 0, WM = 0.0047):
        virtualbata = BetaBranch(self.Zavg, self.Aavg, frac=contribute, I=0, E0=E0, sigma_E0=0, forbiddeness=forbiddeness, WM=WM)
        return virtualbata.BetaSpectrum(x, nu_spectrum=True)

    # function that fit the reference beta spectrum with virtual brances
    def FitData(self, betadata, slicesize):
        self.contribute = {}
        self.E0 = {}
        self.Zlist = {}
        self.Alist = {}
        # fill the sub lists as cached slices
        subx = []
        suby = []
        subyerr = []
        xhigh = betadata.x[-1]
        datacache = np.copy(betadata.y) # preserve the data
        for it in range(len(betadata.x)-1, 0, -1):
            x = betadata.x[it]

            if x < xhigh - slicesize:
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

                    popt, pcov = curve_fit(self.BetaSpectrum, subx, suby, p0 = init_guess, sigma=subyerr, absolute_sigma=True, bounds=(0, [xhigh*1.5, limit]))
                    self.contribute[xhigh] = popt[1]
                    self.E0[xhigh] = popt[0]

                    # subtract the best fit spectrum from beta data
                    for i in range(len(betadata.x)):
                        datacache[i] -= self.BetaSpectrum(betadata.x[i], popt[0], popt[1])

                    #TODO uncertainty process

                    # reset cached slices
                    subx = []
                    suby = []
                    subyerr = []

                xhigh -= slicesize
            else:
                subx.append(x)
                suby.append(datacache[it])
                subyerr.append(betadata.yerr[it])

    # function to calculate summed spectra of virtual branches
    def SumBranches(self, x, thresh = 0, nu_spectrum = True):
        result = 0
        for s in self.E0:
            if s > thresh: # if thresh > 0, look at spectra in selected region
                vb = BetaBranch(self.Zlist[s], self.Alist[s], frac=self.contribute[s], I=0, E0=self.E0[s], sigma_E0=0, forbiddeness=0, WM=0.0047)
                result += vb.BetaSpectrum(x, nu_spectrum = nu_spectrum)
        return result

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
        for istp in self.betadata:
            # define the virtual branches to be fit
            vbnew = VirtualBranch(self.fisIstp[istp])
            vbnew.FitData(self.betadata[istp], slicesize)
            self.vblist[istp] = vbnew

    # Draw plots to test output.
    def DrawVB(self, name, figname, summing = True, setlog = True):
        print("Drawing spectrum...")
        fig = plt.figure()

        xval = np.linspace(0., 10., 200)
        #plt.errorbar(self.betadata[name].x, self.vblist[name].SumBranches(self.betadata[name].x, True))
        for i in range(0, 20):
            plt.errorbar(xval, self.vblist[name].SumBranches(xval, thresh =i*0.5, nu_spectrum = False), fmt='--')
        if setlog: plt.yscale('log')
        plt.errorbar(self.betadata[name].x, self.betadata[name].y, self.betadata[name].yerr, label='beta data')
        plt.errorbar(xval, self.vblist[name].SumBranches(xval, True))

        fig.savefig(figname)


# test
if __name__ == "__main__":
    beta235 = BetaData("./conversionDB/U_235_e_2014.csv")
    #print(beta235.x)
    #print(beta235.y)
    #print(beta235.yerr)

    U235 = FissionIstp(92, 235)
    U235.LoadDB()

    # vbtest = VirtualBranch(U235)
    # vbtest.CalcZAavg(6,7)
    # print(vbtest.Aavg, vbtest.Zavg)

    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)
    convertmodel.VBfit(0.25)
    convertmodel.DrawVB("U235", "U235_convert_test.png")

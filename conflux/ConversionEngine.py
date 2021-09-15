# beta spectra conversion engine:
# input: beta spectra measured on nuclear reactors
# output: neutrino spectra of nuclear reactors

# universal modules
import sys
import csv
import numpy as np
from scipy.optimize import curve_fit

# local modules
from BetaEngine import BetaEngine, BetaBranch
from FPYEngine import FissionModel, FissionIstp

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
    def BetaSpectrum(self, x, contribute, E0, forbiddeness = 0, WM = 0.0047):
        virtualbata = BetaBranch(self.Zavg, self.Aavg, frac=contribute, I=0, E0=E0, sigma_E0=0, forbiddeness=forbiddeness, WM=WM)
        return virtualbata.BetaSpectrum(x)

    # define the theoretical neutrino spectrum shape
    def NueSpectrum(self, x, contribute, E0, forbiddeness = 0, WM = 0.0047):
        virtualbata = BetaBranch(self.Zavg, self.Aavg, frac=contribute, I=0, E0=E0, sigma_E0=0, forbiddeness=forbiddeness, WM=WM)
        return irtualbata.BetaSpectrum(x, nu_spectrum=True)

    def FitData(self, betadata, slicesize):
        self.contribute = {}
        self.E0 = {}
        self.Zlist = {}
        self.Alist = {}
        # fill the sub lists as slices
        subx = []
        suby = []
        subyerr = []
        xhigh = betadata.x[-1]
        for it in range(len(betadata.x)-1, 0, -1):
            x = betadata.x[it]

            if x < xhigh - slicesize:
                # when the sublist is filled in this slice, do fitting
                if len(subx)>1 and len(suby) == len(subx) :
                    self.CalcZAavg(xhigh-slicesize, xhigh)
                    self.Zlist[xhigh] = self.Zavg
                    self.Alist[xhigh] = self.Aavg
                    print(self.Zavg, self.Aavg, subx, suby, subyerr)
                    popt, pcov = curve_fit(self.BetaSpectrum, subx, suby, sigma=subyerr, absolute_sigma=True)
                    self.contribute[xhigh] = popt[0]
                    self.E0[xhigh] = popt[1]
                    print(popt[1])
                    # subtract the best fit spectrum from beta data
                    for i in range(len(betadata.x)):
                        betadata.y[i] - self.BetaSpectrum(popt[0], popt[1])

                    #TODO uncertainty process

                    subx = []
                    suby = []
                    subyerr = []

                xhigh -= slicesize
            else:
                subx.append(x)
                suby.append(betadata.y[it])
                subyerr.append(betadata.yerr[it])

# class that search for best fit vertual branch and calculate total neutrino flux
class ConversionEngine:
    def __init__(self):
        self.betadata = {}
        self.fissionfrac = {}
        self.fisIstp = {}
        self.vblist = {}

    def AddBetaData(self, betadata, fisIstp, name, frac):
        self.betadata[name] = betadata
        self.fissionfrac[name] = frac
        self.fisIstp[name] = fisIstp

    def VBfit(self, slicesize = 0.25):
        for istp in self.betadata:
            # define the virtual branches to be fit
            vbnew = VirtualBranch(self.fisIstp[istp])
            vbnew.FitData(self.betadata[istp], slicesize)
            self.vblist[istp] = vbnew
            print(self.vblist[istp].E0)

# test
if __name__ == "__main__":
    beta235 = BetaData("./conversionDB/U_235_e_2014.csv")
    #print(beta235.x)
    #print(beta235.y)
    #print(beta235.yerr)

    U235 = FissionIstp(92, 235)
    U235.LoadDB()

    vbtest = VirtualBranch(U235)
    vbtest.CalcZAavg(6,7)
    print(vbtest.BetaSpectrum(8, 0.1, 9))
    #print(vbtest.Aavg, vbtest.Zavg)

    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)
    convertmodel.VBfit()

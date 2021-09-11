# beta spectra conversion engine:
# input: beta spectra measured on nuclear reactors
# output: neutrino spectra of nuclear reactors

# universal modules
import sys
import csv
import numpy as np

# local modules
from conflux.BetaEngine import BetaEngine
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
                y = float(row["dNe"])
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
    def CalcZAavg(self, E0, dE0):
        Elow = E0-dE0
        Ehigh = E0+dE0
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

    def BetaSpectrum(self, E0, forbiddeness = 0, WM = 0.0047):
        return 0


class ConversionEngine:
    def __init__():
        betadata

# test
if __name__ == "__main__":
    beta235 = BetaData("./conversionDB/U_235_e_2014.csv")
    print(beta235.x)
    print(beta235.y)
    print(beta235.yerr)

    U235 = FissionIstp(92, 235)
    U235.LoadDB()

    vbtest = VirtualBranch(U235)
    vbtest.CalcZAavg()
    print(vbtest.Aavg, vbtest.Zavg)

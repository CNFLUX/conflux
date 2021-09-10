# beta spectra conversion engine:
# input: beta spectra measured on nuclear reactors
# output: neutrino spectra of nuclear reactors

# universal modules
import sys
import csv

# local modules
from conflux.BetaEngine import BetaEngine

# Class that reads conversion DB and form reference spectra
class betadata:
    def __init__(self, inputDB):
        self.x = []
        self.y = []
        self.yerr = []
        self.inputDB = inputDB
        self.LoadDB(inputDB)

    def LoadDB(self, inputDB):
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

                self.x.append(E)
                self.y.append(float(row["Ne"]))
                self.yerr.append(float(row["dNe"]))


class ConversionEngine:
    def __init__():
        betadata

# test
if __name__ == "__main__":
    beta235 = betadata("./conversionDB/U_235_e_2014.csv")
    print(beta235.x)
    print(beta235.y)
    print(beta235.yerr)

# fission product and spectra summation engine

# universal modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

# local modules
from BetaEngine import BetaEngine
from FPYEngine import FissionModel, FissionIstp

class SumEngine:
    def __init__(self, neutrino=True):
        self.FPYlist = {}
        self.betaSpectraList = {}
        self.neutrino = neutrino # neutrino or electon spectrum

    # method to add fission/non-fissile/non-equilibrium isotopes into the engine
    def AddModel(self, fissionModel, W=1.0):
        for FPZAI in fissionModel.FPYlist:
            if FPZAI not in self.FPYlist:
                self.FPYlist[FPZAI] = fissionModel.FPYlist[FPZAI]
                self.FPYlist[FPZAI].y *= W
                self.FPYlist[FPZAI].yerr *= W
            else:
                self.FPYlist[FPZAI].y + fissionModel.FPYlist[FPZAI].y*W
                self.FPYlist[FPZAI].yerr + fissionModel.FPYlist[FPZAI].yerr*W

    def NormalizeFP(self):
        self.sum = 0
        for FPZAI in self.FPYlist:
            self.sum += self.FPYlist[FPZAI].y
        for FPZAI in self.FPYlist:
            self.FPYlist[FPZAI].y /= self.sum
            self.FPYlist[FPZAI].yerr /= self.sum

    def CalcReactorSpectrum(self, betaSpectraDB, DB='betaDB/betaDB.xml', binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0):
        bins = int(erange/binwidths)
        self.reactorSpectrum = np.zeros(bins)
        self.spectrumErr = np.zeros(bins)
        self.bins = np.arange(0, erange, binwidths)
        self.missingBranch = []
        self.missingCont = 0.0

        for FPZAI in self.FPYlist:
            if FPZAI in betaSpectraDB.spectralist:
                self.betaSpectraList[FPZAI] = betaSpectraDB.spectralist[FPZAI]
                self.reactorSpectrum += self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y
                self.spectrumErr += self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].yerr

            # save the list of missing branches
            else:
                self.missingCont += self.FPYlist[FPZAI].y
                self.missingBranch.append(FPZAI)

    def Draw(self, figname, summing = True, logy=True, frac=False):
        print("Drawing spectrum...")
        fig, ax = plt.subplots()
        if (not logy):
            ax.set_xlim([0, 10])
            ax.set_ylim([0, 1])
            if (summing == True):
                ax.plot(self.bins, self.reactorSpectrum)
            if (frac == True):
                lines = []
                labels = []
                for FPZAI in self.betaSpectraList:
                    lines += ax.plot(self.bins, self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y)
                    labels.append(str(FPZAI))
                    # if (self.betaSpectraList[FPZAI][60]>0 and np.max(self.betaSpectraList[FPZAI])> 5e-3):
                    #     print(FPZAI, self.betaSpectraList[FPZAI][60], self.FPYlist[FPZAI].y)
                plt.legend(lines, labels)
                ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='neutrino spectrum')
        else:
            ax.set_xlim([0, 15])
            ax.set_ylim([1e-7, 10])
            if (summing == True):
                ax.semilogy(self.bins, self.reactorSpectrum)

            if (frac == True):
                for FPZAI in self.betaSpectraList:
                    ax.semilogy(self.bins, self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y)
                    ax.legend(str(FPZAI))

                    # if (self.betaSpectraList[FPZAI][60]>0 and np.max(self.betaSpectraList[FPZAI])> 5e-3):
                    #     print(FPZAI, self.betaSpectraList[FPZAI][60], self.FPYlist[FPZAI].y)
            ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='neutrino spectrum')
        fig.savefig(figname)

if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadDB()
    U238 = FissionIstp(92, 238)
    U238.LoadDB()
    Pu239 = FissionIstp(94, 239)
    Pu239.LoadDB()

    model = FissionModel()
    #
    model.AddContribution(isotope=U235, Ei = 0, fraction=1)
    #model.AddContribution(isotope=U238, Ei = 0.5, fraction=0.5)
    #model.AddContribution(isotope=Pu239, Ei = 0, fraction=0.5)
    #model.AddIstp(39, 96, 1.0)

    result = SumEngine()
    result.AddModel(model)

    betaSpectraDB = BetaEngine(result.FPYlist.keys())
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0)

    result.CalcReactorSpectrum(betaSpectraDB)
    print(result.reactorSpectrum)
    result.Draw("test_spectrum_test.png")
    with open("U235test.csv", "w") as output:
        write = csv.writer(output)
        write.writerow(result.reactorSpectrum)
    #print(result.missingCont)
    #print(result.missingBranch)

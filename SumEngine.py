# fission product and spectra summation engine

# universal modules
import sys
import numpy as np
import matplotlib.pyplot as plt

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

    def CalcReactorSpectrum(self, DB='betaDB/betaDB.xml', binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0):
        bins = int(erange/binwidths)
        self.reactorSpectrum = np.zeros(bins)
        self.bins = np.arange(0, erange, binwidths)

        self.betaSpectraDB = BetaEngine(self.FPYlist.keys())
        self.betaSpectraDB.CalcBetaSpectra(DB, nu_spectrum=self.neutrino, binwidths=binwidths, lower=lower, thresh=thresh, erange = erange)
        for FPZAI in self.FPYlist:
            if FPZAI in self.betaSpectraDB.spectralist:
                #print(FPZAI, self.betaSpectraDB.spectralist[FPZAI])
                self.betaSpectraList[FPZAI] = self.betaSpectraDB.spectralist[FPZAI]
                self.reactorSpectrum += self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y

    def Draw(self, figname, logy=True):
        print("Drawing spectrum...")
        fig, ax = plt.subplots()
        if (not logy):
            ax.plot(self.bins, self.reactorSpectrum)
        else:
            ax.semilogy(self.bins, self.reactorSpectrum)
        ax.set(xlabel='E (MeV)', ylabel='I', title='reactor neutrino spectrum')
        fig.savefig(figname)



if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadDB()
    Pu239 = FissionIstp(94, 239)
    Pu239.LoadDB()

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0, fraction=0.6, d_frac=0.05)
    model.AddContribution(isotope=Pu239, Ei = 0, fraction=0.4, d_frac=0.05)

    result = SumEngine()
    result.AddModel(model)
    result.CalcReactorSpectrum()
    print(result.reactorSpectrum)
    result.Draw("reactor_spectrum_test.png")

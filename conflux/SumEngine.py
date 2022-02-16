# fission product and spectra summation engine

# universal modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

# local modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp

class SumEngine:
    def __init__(self, neutrino=True):
        self.FPYlist = {}
        self.betaSpectraList = {}
        self.betaUncertainty = {}
        self.neutrino = neutrino # neutrino or electon spectrum

    def Clear(self):
        self.FPYlist = {}
        self.betaUncertainty = {}
        self.betaSpectraList = {}

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
                self.FPYlist[FPZAI].AddCovariance(fissionModel.FPYlist[FPZAI])

    def NormalizeFP(self):
        self.sum = 0
        for FPZAI in self.FPYlist:
            self.sum += self.FPYlist[FPZAI].y
        for FPZAI in self.FPYlist:
            self.FPYlist[FPZAI].y /= self.sum
            self.FPYlist[FPZAI].yerr /= self.sum

    def CalcReactorSpectrum(self, betaSpectraDB, binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0):
        print("calculating beta spectra...")
        bins = int(erange/binwidths)
        self.reactorSpectrum = np.zeros(bins)
        self.spectrumUnc = np.zeros(bins)   # total uncertainty
        self.modelUnc = np.zeros(bins)
        self.yieldUnc = np.zeros(bins)
        self.bins = np.arange(0, erange, binwidths)
        self.missingBranch = []
        self.missingCont = 0.0

        for FPZAI in self.FPYlist:
            if FPZAI in betaSpectraDB.spectralist:
                self.betaSpectraList[FPZAI] = betaSpectraDB.spectralist[FPZAI]
                self.betaUncertainty[FPZAI] = betaSpectraDB.uncertaintylist[FPZAI]
                self.reactorSpectrum += self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y
                self.yieldUnc += self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].yerr
                self.modelUnc += self.betaUncertainty[FPZAI]*self.FPYlist[FPZAI].y
                #self.spectrumUnc = self.yieldUnc+self.modelUnc

            # save the list of missing branches
            else:
                self.missingCont += self.FPYlist[FPZAI].y
                self.missingBranch.append(FPZAI)

        # Uncertainty calculation
        for i in self.FPYlist:
            sigmai = 0
            for j in self.FPYlist:
                if i in betaSpectraDB.spectralist and j in betaSpectraDB.spectralist:
                    yi = self.FPYlist[i].y
                    yerri = self.FPYlist[i].yerr
                    fi = self.betaSpectraList[i]

                    yj = self.FPYlist[j].y
                    yerrj = self.FPYlist[j].yerr
                    fj = self.betaSpectraList[j]

                    sigmay_ij = self.FPYlist[i].cov[j] # FIXME: Why is sigmay_ij different from desired covariance elements? Where is the
                    #sigmai += self.FPYlist[i].cov[j]
                    #if (i == j): print(i, sigmay_ij, self.FPYlist[i].cov[i])
                    self.spectrumUnc += fi*sigmay_ij*fj

            #print("cov sum", i, sigmai)

        self.spectrumUnc = np.sqrt(self.spectrumUnc)

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
                    if (self.betaSpectraList[FPZAI][90]>0):
                        lines += ax.plot(self.bins, self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y)
                        labels.append(str(FPZAI))
                plt.legend(lines, labels)
                ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='neutrino spectrum')
        else:
            #ax.set_xlim([0, 15])
            ax.set_xlim([1e-7, 10])
            ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='neutrino spectrum')

            if (frac == True):
                lines = []
                labels = []
                for FPZAI in self.betaSpectraList:
                    if (self.betaSpectraList[FPZAI][10]>0):
                        #print(FPZAI, self.FPYlist[FPZAI].y, self.betaSpectraList[FPZAI][90])
                        lines += ax.plot(self.bins, \
                        self.betaSpectraList[FPZAI]*self.FPYlist[FPZAI].y)
                        labels.append(str(FPZAI))

            if (summing == True):
                ax.semilogy(self.bins, self.reactorSpectrum)

        fig.savefig(figname)

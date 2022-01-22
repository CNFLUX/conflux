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
                print("FPYLIST_COV:", FPZAI, self.FPYlist[FPZAI].cov[FPZAI])
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

    def CalcReactorSpectrum(self, betaSpectraDB, DB='betaDB/betaDB.xml', binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0):
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

                    sigmay_ij = yerri*yerrj*self.FPYlist[i].cov[j] # FIXME: Why is sigmay_ij different from desired covariance elements? Where is the
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

if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadDB()
    U235.LoadCovariance()
    # U238 = FissionIstp(92, 238)
    # U238.LoadDB()
    # U233 = FissionIstp(92, 233)
    # U233.LoadDB()
    # U234 = FissionIstp(92, 234)
    # U234.LoadDB()
    # Pu239 = FissionIstp(94, 239)
    # Pu239.LoadDB()
    # Pu241 = FissionIstp(94, 241)
    # Pu241.LoadDB()

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0, fraction=1)
    #model.AddContribution(isotope=U234, Ei = 0.5, fraction=1)
    #model.AddContribution(isotope=U233, Ei = 0, fraction=1)
    #model.AddContribution(isotope=Pu241, Ei = 0, fraction=0.0572)
    #model.AddIstp(39, 96, 1.0)

    result = SumEngine()
    result.AddModel(model)

    betaSpectraDB = BetaEngine(result.FPYlist.keys())
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1, lower=-1.0, thresh=0.0, erange = 10.0)

    result.CalcReactorSpectrum(betaSpectraDB, erange = 10.0)
    summed_spect = result.reactorSpectrum
    summed_err = result.spectrumUnc
    summed_model_err = result.modelUnc

    result.Draw("Commercial.png", frac=False)
    result.Clear()

    # model1 = FissionModel()
    # model1.AddContribution(isotope=U235, Ei = 0, fraction=1)
    # result.AddModel(model1)
    # result.CalcReactorSpectrum(betaSpectraDB, erange = 10.0)
    # m1_spect = result.reactorSpectrum
    # result.Clear()
    #
    # model2 = FissionModel()
    # model2.AddContribution(isotope=U234, Ei = 0.5, fraction=1)
    # result.AddModel(model2)
    # result.CalcReactorSpectrum(betaSpectraDB, erange = 10.0)
    # m2_spect = result.reactorSpectrum
    # result.Clear()
    #
    # model3 = FissionModel()
    # model3.AddContribution(isotope=U233, Ei = 0, fraction=1)
    # result.AddModel(model3)
    # result.CalcReactorSpectrum(betaSpectraDB, erange = 10.0)
    # m3_spect = result.reactorSpectrum
    # result.Clear()
    #
    # model4 = FissionModel()
    # model4.AddContribution(isotope=U235, Ei = 0, fraction=1)
    # result.AddModel(model4)
    # result.CalcReactorSpectrum(betaSpectraDB, erange = 10.0)
    # m4_spect = result.reactorSpectrum
    # result.Clear()

    fig, ax = plt.subplots()
    ax.set_ylim([1e-6, 10])
    ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    ax.fill_between(result.bins, summed_spect+summed_err, summed_spect-summed_err, alpha=.5, linewidth=0)
    ax.plot(result.bins, summed_spect, label="Summed")
    # ax.errorbar(result.bins, summed_spect, yerr = summed_model_err, label="Beta model uncertainty")
    ax.legend()
    #ax.semilogy(result.bins, m4_spect, label="U235")
    #ax.legend()
    # ax.semilogy(result.bins, m2_spect, label="U234")
    # ax.legend()
    # ax.semilogy(result.bins, m3_spect, label="U233")
    # ax.legend()
    # ax.semilogy(result.bins, m4_spect, label="Pu241")
    # ax.legend()
    fig.savefig("drawtest1.png")

    with open("Commercial.csv", "w") as output:
        write = csv.writer(output)
        #write.writerow(result.reactorSpectrum)
    #print(result.missingCont)
    #print(result.missingBranch)

import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import operator

# conflux modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine

if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB(defaultDB='ENDF')
    U235.LoadCorrelation(defaultDB='ENDF')

        # Pu239 = FissionIstp(94, 239)
        # Pu239.LoadFissionDB()
        # Pu239.LoadCorrelation()
    #U235.CalcCovariance(Ei=0)

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0, fraction=1)
    #model.AddContribution(isotope=Pu239, Ei = 0, fraction=0)
    #model.AddContribution(isotope=U233, Ei = 0, fraction=1)
    #model.AddContribution(isotope=Pu241, Ei = 0, fraction=0.0572)
    #model.AddIstp(39, 96, 1.0)

    result = SumEngine()
    result.AddModel(model)

    betaSpectraDB = BetaEngine(result.FPYlist.keys(), binwidths=0.1, spectRange=[0.0, 15.0])
    #betaSpectraDB = BetaEngine(newlist)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])
        
    betaDBBase = BetaEngine()
    count = 0
    newlist = []
    for nuclide in (sorted(U235.CFPY[0].values(), key=operator.attrgetter('y'), reverse=True)):
        if nuclide.FPZAI in betaDBBase.istplist:
            fig, ax = plt.subplots()
            ax.set(xlabel='E (MeV)', ylabel='variance-covariance', title='neutrino spectrum uncertainty')
            yerr = nuclide.yerr
            ax.plot(betaSpectraDB.xbins, yerr**2*betaSpectraDB.istplist[nuclide.FPZAI].spectrum, label=betaSpectraDB.istplist[nuclide.FPZAI].name+" variance")
            positive=np.zeros(len(betaSpectraDB.xbins))
            negative=np.zeros(len(betaSpectraDB.xbins))
            for FPZAI, frac in (sorted(nuclide.cov.items(), key=lambda item: abs(item[1]), reverse=True)):
                if FPZAI in betaSpectraDB.istplist:
                    if frac >0:
                        positive += frac*betaSpectraDB.istplist[FPZAI].spectrum
                        ax.fill_between(betaSpectraDB.xbins, positive,  positive-frac*betaSpectraDB.istplist[FPZAI].spectrum, alpha = 0.4)
                    else:
                        negative += frac*betaSpectraDB.istplist[FPZAI].spectrum
                        ax.fill_between(betaSpectraDB.xbins, negative,  negative-frac*betaSpectraDB.istplist[FPZAI].spectrum, alpha = 0.4)
            ax.plot(betaSpectraDB.xbins, positive+negative, label=betaSpectraDB.istplist[nuclide.FPZAI].name+" covariance")
            ax.legend()
            fig.savefig(betaDBBase.istplist[nuclide.FPZAI].name+"_cov.png")
            newlist.append(nuclide.FPZAI)
            count += 1
        if count == 5:
            break

    result.CalcReactorSpectrum(betaSpectraDB, spectRange=[0.0, 15.0], branchErange=[0.0, 20.0], processMissing=False)
    summed_spect = result.reactorSpectrum
    summed_err = result.spectrumUnc
    summed_model_err = result.modelUnc
    summed_yerr = result.yieldUnc
    
    print(result.totalYield)
    print(result.missingCount)
    print(result.missingBranch)
    
    result.CalcReactorSpectrum(betaSpectraDB, spectRange=[0.0, 15.0], branchErange=[0.0, 20.0], processMissing=True)
    miss_spect = result.reactorSpectrum
    miss_err = result.spectrumUnc
    miss_model_err = result.modelUnc
    miss_yerr = result.yieldUnc
    
    # print(result.totalYield)
    # print(result.missingCount)
    # print(result.missingBranch)
    
    result.Clear()

    fig, ax = plt.subplots()
    #ax.set_ylim([-1, 1])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    ax.fill_between(result.bins, miss_spect+miss_yerr, miss_spect-miss_yerr, alpha=.5, linewidth=0, label="fission product error")
    ax.fill_between(result.bins, miss_spect+miss_model_err, miss_spect-miss_model_err, alpha=.5, linewidth=0, label="beta model error")
    ax.plot(result.bins, miss_spect, label="w/ miss info")
    ax.plot(result.bins, miss_spect-summed_spect, label="missing info")
    # ax.plot(result.bins, miss_spect, label="w/ miss info")
    # ax.fill_between(result.bins, summed_err, -summed_err, alpha=.5, linewidth=0)
    # ax.fill_between(result.bins, summed_yerr, -summed_yerr, alpha=.5, linewidth=0)
    # ax.errorbar(result.bins, summed_spect, yerr = summed_model_err, label="Beta model uncertainty")
    ax.legend()

    fig.savefig("235U_ENDF_TOP_linear.png")
    

    fig, ax = plt.subplots()
    ax.set_xlim([0, 10])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', ylabel='relative error (%)')
    ax.fill_between(result.bins, miss_yerr/miss_spect*100, -miss_yerr/miss_spect*100, label="fission product error", alpha=.5)
    ax.fill_between(result.bins, miss_model_err/miss_spect*100, -miss_model_err/miss_spect*100, label="beta model error", alpha=.5)
    ax.plot(result.bins, miss_spect-miss_spect)
    #ax.plot(result.bins, miss_spect-summed_spect, label="missing info")
    # ax.plot(result.bins, miss_spect, label="w/ miss info")
    # ax.fill_between(result.bins, summed_err, -summed_err, alpha=.5, linewidth=0)
    # ax.fill_between(result.bins, summed_yerr, -summed_yerr, alpha=.5, linewidth=0)
    # ax.errorbar(result.bins, summed_spect, yerr = summed_model_err, label="Beta model uncertainty")
    ax.legend()

    fig.savefig("235U_ENDF_Unc.png")
    
    fig, ax = plt.subplots()
    #ax.set_ylim([-1, 1])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', ylabel='delta neutrino/decay/MeV', title='U-235 neutrino flux')
    ax.plot(result.bins, miss_spect-summed_spect)

    ax.legend()

    fig.savefig("235_239_Missing.png")


    with open("Commercial.csv", "w") as output:
        write = csv.writer(output)

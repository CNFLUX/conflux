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
    xbins = np.arange(0, 20, 0.1)

    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB(defaultDB='JEFF')
    #U235.LoadCorrelation(defaultDB='ENDF')

    Pu239 = FissionIstp(94, 239)
    Pu239.LoadFissionDB(defaultDB='JEFF')
    # Pu239.LoadCorrelation()
    #U235.CalcCovariance(Ei=0)

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 14, fraction=1, IFP=True)
    model.SaveToFile('FPY_235_JEFF_IFP_14MeV.csv')
    # model.AddContribution(isotope=Pu239, Ei = 0.4, fraction=1)
    # model.SaveToFile('FPY_239_JEFF_IFP.csv')

    #model.AddContribution(isotope=U233, Ei = 0, fraction=1)
    #model.AddContribution(isotope=Pu241, Ei = 0, fraction=0.0572)
    #model.AddIstp(39, 96, 1.0)

    sum1 = SumEngine(xbins = xbins)
    sum1.AddModel(model)

    betaSpectraDB = BetaEngine(sum1.FPYlist.keys(), xbins=xbins)
    #betaSpectraDB = BetaEngine(newlist)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    sum1.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False,  ifp_begin = 10,  ifp_end = 100)
    summed_spect = sum1.spectrum
    summed_err = sum1.uncertainty
    summed_model_err = sum1.modelUnc
    summed_yerr = sum1.yieldUnc

    #result.Clear()

    sum2 = SumEngine(xbins=xbins)
    sum2.AddModel(model)
    sum2.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False,  ifp_begin = 0,  ifp_end = 10)
    miss_spect = sum2.spectrum
    miss_err = sum2.uncertainty
    miss_model_err = sum2.modelUnc
    miss_yerr = sum2.yieldUnc
    sum2.SaveToFile('235U_nu_jeff_ifp.csv')


    sum2.Clear()

    fig, ax = plt.subplots()
    #ax.set_ylim([-1, 1])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    # ax.fill_between(sum2.xbins, miss_spect+miss_yerr, miss_spect-miss_yerr, alpha=.5, linewidth=0, label="fission product error")
    # ax.fill_between(sum2.xbins, miss_spect+miss_model_err, miss_spect-miss_model_err, alpha=.5, linewidth=0, label="beta model error")
    ax.plot(sum2.xbins, miss_spect, label="10-100")
    #ax.plot(sum2.xbins, summed_spect, label="w/o info")

    ax.plot(sum2.xbins, summed_spect, label="0-10")

    ax.legend()

    fig.savefig("235U_ENDF_jeff_14_MeV.png")


    fig, ax = plt.subplots()
    ax.set_xlim([0, 10])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', ylabel='relative error (%)')
    ax.fill_between(sum2.xbins, miss_yerr/miss_spect*100, -miss_yerr/miss_spect*100, label="fission product error", alpha=.5)
    ax.fill_between(sum2.xbins, miss_model_err/miss_spect*100, -miss_model_err/miss_spect*100, label="beta model error", alpha=.5)
    ax.plot(sum2.xbins, miss_spect-miss_spect)
    #ax.plot(result.bins, miss_spect-summed_spect, label="missing info")
    # ax.plot(result.bins, miss_spect, label="w/ miss info")
    # ax.fill_between(result.bins, summed_err, -summed_err, alpha=.5, linewidth=0)
    # ax.fill_between(result.bins, summed_yerr, -summed_yerr, alpha=.5, linewidth=0)
    # ax.errorbar(result.bins, summed_spect, yerr = summed_model_err, label="Beta model uncertainty")
    ax.legend()

    fig.savefig("235U_ENDF_jeff_14_MeV_unc.png")

    fig, ax = plt.subplots()
    #ax.set_ylim([-1, 1])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', ylabel='delta neutrino/decay/MeV', title='U-235 neutrino flux')
    ax.plot(sum2.xbins, miss_spect-summed_spect)

    ax.legend()

    fig.savefig("235_239_Missing_jeff.png")


    with open("Commercial.csv", "w") as output:
        write = csv.writer(output)

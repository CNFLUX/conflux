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

    U235 = FissionIstp(92, 235, Ei = 0.4, DB='JEFF', IFPY=True)
    U235.LoadFissionDB(Ei = 0.5)
    U235.LoadCorrelation(DB='JEFF')

    Pu239 = FissionIstp(94, 239, Ei = 0)
    Pu239.LoadFissionDB()
    Pu239.LoadCorrelation()
    U235.CalcCovariance()

    model = FissionModel()
    model.AddContribution(isotope=U235, fraction=1)
    model.SaveToFile('FPY_235_JEFF.csv')
    #model.AddContribution(isotope=Pu239, Ei = 0, fraction=0)
    #model.AddContribution(isotope=U233, Ei = 0, fraction=1)
    #model.AddContribution(isotope=Pu241, Ei = 0, fraction=0.0572)
    #model.AddIstp(39, 96, 1.0)

    sum1 = SumEngine(xbins = xbins)
    sum1.AddModel(model)

    betaSpectraDB = BetaEngine(xbins=xbins)
    #betaSpectraDB = BetaEngine(newlist)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    U235.CalcBetaSpectra(betaSpectraDB)
    print(U235.spectrum)
    betaDBBase = BetaEngine()
    count = 0
    newlist = []
    for nuclide in (sorted(U235.FPYlist.values(), key=operator.attrgetter('y'), reverse=True)):
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

    sum1.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False)
    summed_spect = sum1.spectrum
    summed_err = sum1.uncertainty
    summed_model_err = sum1.modelUnc
    summed_yerr = sum1.yieldUnc

    sum2 = SumEngine(xbins=xbins)
    sum2.AddModel(model)
    sum2.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=True)
    miss_spect = sum2.spectrum
    miss_err = sum2.uncertainty
    miss_model_err = sum2.modelUnc
    miss_yerr = sum2.yieldUnc
    sum2.SaveToFile('235U_nu_endf_test.csv')
    sum1.SaveToFile('235U_nu_jeff_nomiss_test.csv')

    sum2.Clear()

    fig, ax = plt.subplots()
    ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    ax.fill_between(sum2.xbins, miss_spect+miss_yerr, miss_spect-miss_yerr, alpha=.5, linewidth=0, label="fission product error")
    ax.fill_between(sum2.xbins, miss_spect+miss_model_err, miss_spect-miss_model_err, alpha=.5, linewidth=0, label="beta model error")
    ax.plot(sum2.xbins, miss_spect, label="w/ miss info")

    ax.plot(sum2.xbins, miss_spect-summed_spect, label="missing info")
    ax.errorbar(sum2.xbins, summed_spect, yerr = summed_model_err, label="Beta model uncertainty")
    ax.legend()

    fig.savefig("235U_ENDF_TOP_linear_jeff_test.png")

    fig, ax = plt.subplots()
    ax.set_xlim([0, 10])
    ax.set(xlabel='E (MeV)', ylabel='relative error (%)')
    ax.fill_between(sum2.xbins, miss_yerr/miss_spect*100, -miss_yerr/miss_spect*100, label="fission product error", alpha=.5)
    ax.fill_between(sum2.xbins, miss_model_err/miss_spect*100, -miss_model_err/miss_spect*100, label="beta model error", alpha=.5)
    ax.plot(sum2.xbins, miss_spect-miss_spect)
    ax.legend()

    fig.savefig("235U_ENDF_Unc_jeff_test.png")

    fig, ax = plt.subplots()
    ax.set(xlabel='E (MeV)', ylabel='delta neutrino/decay/MeV', title='U-235 neutrino flux')
    ax.plot(sum2.xbins, miss_spect-summed_spect)
    ax.legend()
    fig.savefig("235_239_Missing_jeff_test.png")

    with open("Commercial.csv", "w") as output:
        write = csv.writer(output)

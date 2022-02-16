import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

# conflux modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine

if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB()
    U235.LoadCorrelation()
    U235.CalcCovariance(Ei=0)

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
    summed_yerr = result.yieldUnc

    result.Draw("Commercial.png", frac=False)
    result.Clear()


    fig, ax = plt.subplots()
    ax.set_ylim([-1, 1])
    ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    #ax.fill_between(result.bins, summed_spect+summed_err, summed_spect-summed_err, alpha=.5, linewidth=0)
    #ax.plot(result.bins, summed_spect, label="Summed")
    ax.fill_between(result.bins, summed_err, -summed_err, alpha=.5, linewidth=0)
    ax.fill_between(result.bins, summed_yerr, -summed_yerr, alpha=.5, linewidth=0)
    # ax.errorbar(result.bins, summed_spect, yerr = summed_model_err, label="Beta model uncertainty")
    ax.legend()

    fig.savefig("uncertainty.png")

    with open("Commercial.csv", "w") as output:
        write = csv.writer(output)

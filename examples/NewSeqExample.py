import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import operator

# conflux modules
from conflux.BetaEngine import BetaEngine, BetaIstp
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine

if __name__ == "__main__":
    # Define spectrum binning
    xbins = np.arange(0, 10, 0.1)

    # Calculate individual beta spect
    betaSpectraDB = BetaEngine(xbins=xbins)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    Pu239 = FissionIstp(94, 239, Ei = 0)
    Pu239.LoadFissionDB()
    Pu239.LoadCorrelation()
    Pu239.CalcBetaSpectra(betaSpectraDB)

    # Define and calculate fission and nu spect
    U235 = FissionIstp(92, 235, Ei = 0.5, DB='ENDF', IFPY=False)
    U235.LoadFissionDB(Ei = 0.5)
    U235.LoadCorrelation(DB='ENDF')
    U235.CalcBetaSpectra(betaSpectraDB)

    # Define and calculate fission and nu spect
    U235T = FissionIstp(92, 235, Ei = 0, DB='ENDF', IFPY=False)
    U235T.LoadFissionDB()
    U235T.LoadCorrelation(DB='ENDF')
    U235T.CalcBetaSpectra(betaSpectraDB)



    # Define summation engine
    newsum = SumEngine(betaSpectraDB)

    # Add fission spect
    newsum.AddFissionIstp(U235, "U235", 100, 0)
    # newsum.AddFissionIstp(Pu239, "Pu239", 50, 7)

    pu241beta = BetaIstp(94, 241, 0, 21.5, 200, "Pu241beta", xbins=xbins)

    # Add beta spect
    # newsum.AddBetaIstp(betaSpectraDB.istplist[922390], "U239", count = 750, d_count = 2)
    # newsum.AddBetaIstp(pu241beta, "Pu241beta", count = 40, d_count = 2)

    # sum all spect
    newsum.CalcReactorSpectrum()

    # final output
    finaloutput = newsum.spectrum

    # # edit contributions freely
    # for i in range(50):
    #     u239number = 750-2*i
    #     newsum.EditContribution("U239", count = u239number, d_count = 2)
    #     newsum.EditContribution("Pu241beta", count = 40-i, d_count = 2)
    #     # recalc spect each time
    #     newsum.CalcReactorSpectrum()

    fig, ax = plt.subplots()
    ax.set_xlim([0, 10])
    ax.set(xlabel='E (MeV)', ylabel='neutrino count')
    ax.errorbar(newsum.xbins, newsum.spectrum, yerr=newsum.uncertainty, label="test spectrum")
    ax.legend()
    plt.show()

from conflux.FPYEngine import FissionIstp, FissionModel
from conflux.ConversionEngine import ConversionEngine, BetaData, VirtualBranch
from conflux.SumEngine import SumEngine
from conflux.BetaEngine import BetaEngine
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    #This is an example to calculate the spectrum of a unique reactor type
    #Using only the Summation Engine. To see how the summation and conversion
    #Engine can be used side-by-side, please see HybridModelExample.py

    #Initialize FissionIstps U235, Pu239, and Pu241 with an incident
    # neutron energy of 0. Load their databases and their covariances as well.

    U235 = FissionIstp(92,235, 0)
    U235.LoadFissionDB()
    U235.LoadCorrelation()

    Pu239 = FissionIstp(94,239, 0)
    Pu239.LoadFissionDB()
    Pu239.LoadCorrelation()

    Pu241 = FissionIstp(94,241, 0)
    Pu241.LoadFissionDB()
    Pu241.LoadCorrelation()

    #I will also go ahead and define my energy range to be from 0 to 15 MeV with
    #100 keV bins

    e = np.arange(0.,15.,0.1)

    #I will go ahead and load up a BetaEngine to Calculate my beta shapes, and
    #Then I will Calculate the BetaSpectra of each of my individual fissionIstps

    BetaEngineDB = BetaEngine(xbins=e)
    BetaEngineDB.CalcBetaSpectra(branchErange=[0., 15])

    U235.CalcBetaSpectra(BetaEngineDB)
    Pu239.CalcBetaSpectra(BetaEngineDB)
    Pu241.CalcBetaSpectra(BetaEngineDB)

    #Next, I will initialize a Summation engine and add all my Fission Isotopes to it

    SummationEngine = SumEngine(BetaEngineDB)
    SummationEngine.AddFissionIstp(U235, "U235", count = 1)
    SummationEngine.AddFissionIstp(Pu239, "Pu239", count = 1)
    SummationEngine.AddFissionIstp(Pu241, "Pu241", count = 1)

    #Finally, I will go ahead and calculate the total reactor spectrum. Note, that because
    #Of my count number, each of these fission isotopes will contribute the same amount
    #To the total spectrum

    SummationEngine.CalcReactorSpectrum()

    #Finally, I want to plot all of this, so I will do that below.

    Spectrum = SummationEngine.spectrum
    Uncertainty = SummationEngine.uncertainty

    fig = plt.plot()
    plt.errorbar(e, Spectrum, Uncertainty, label="Total Spectrum")
    plt.errorbar(e, U235.spectrum, U235.uncertainty, label = "U235 spec", fmt="--")
    plt.errorbar(e, Pu239.spectrum, Pu239.uncertainty, label = "Pu239 spec", fmt="--")
    plt.errorbar(e, Pu241.spectrum, Pu241.uncertainty, label = "Pu241 spec", fmt="--")


    plt.yscale("log")
    plt.xlabel("Energy (in MeV)")
    plt.ylabel("Spectrum (neutrinos/MeV/Fission)")
    plt.legend()
    plt.savefig("NovelReactor.png")
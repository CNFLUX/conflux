from conflux.BetaEngine import BetaEngine
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # This is the example to run the Beta Engine. The Beta Engine
    # Holds all the spectral shape data of our reaction. As such,
    # I only need to calculate it once (Twice, if I'm switching between
    # Neutrino and Beta spectra). 

    # I'll also go ahead and define an energy range for my calculations
    # In this case, it's from 0 to 8 MeV with 100 keV bins
    e = np.arange(0, 8., 0.1)

    # I'll also start by initializing my BetaEngine with the given energy range
    BetaEngineDB = BetaEngine(xbins = e)


    # Once I've loaded a BetaEngine, I can go ahead and load the beta data
    # Needed to do these calculations. I can adjust the database I use for
    # the nuclear information by specifying a :targetDB" that is not default
    # When loading the Beta database.
    BetaEngineDB.LoadBetaDB()

    # Finally, I can calculate the beta spectrum of all isotopes inside my 
    # Target Database. I specify a branchErange, which adjusts what isotopes
    # my engine will calculate shapes for (isotopes whose branches have endpoint
    # energies that fall outside the range will not have their beta shapes calculated)
    # Also note, I can again calculate the spectra for neutrinos as well as betas
    BetaEngineDB.CalcBetaSpectra(branchErange=[0,8], nu_spectrum=True)

    # Each BetaEngine has a dictionary of beta isotopes, which is where the spectral information
    # is stored. That means, I can iterate through the list of Beta Isotopes, and pick one
    # That I want to plot, or look at more spectral information of. Again, each isotope
    # has a unique FPZAI number. See FissionIstp_Example.py for more information.
    for FPZAI in BetaEngineDB.istplist:
        print(BetaEngineDB.istplist[FPZAI].ZAI, BetaEngineDB.istplist[FPZAI].name)


    # I will also go ahead and plot the spectra and uncertainties of Cs-141 using its'
    # Unique FPZAI number.
    cs141spectrum = BetaEngineDB.istplist[551410].spectrum
    cs141uncertainty = BetaEngineDB.istplist[551410].uncertainty
    plt.figure()
    plt.errorbar(BetaEngineDB.xbins, cs141spectrum, cs141uncertainty, label="Cs-141 Spectrum")
    plt.legend()
    plt.xlabel("E (MeV)")
    plt.ylabel("neutrino/MeV")
    plt.savefig("beta_spectrum.pdf")

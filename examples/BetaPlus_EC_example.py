import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaPlusEngine import BetaPlusEngine

if __name__ == "__main__":
    #One thing I didn't specify in the BetaEngine_Example is that you can in fact input a 
    #Custom list of isotopes whose Beta Spectra you want to calculate. You do not need to 
    #Load and Calculate every single Fission Product, and for certain calculations, it's
    #Less resource intensive to input a list.

    #The reaction starts in much the same way as it did in the BetaEngine_Example.py file, 
    #Initialize some energy range, load up a BetaEngine, and then calculate the BetaSpectra
    #However, before I Initialize the BetaEngine, I will define a list of isotopes by their
    #ZAI number. These isotopes will be the ones whose betaSpectrum I want to calculate.

    #initialize energy
    e = np.arange(0.,3.,0.01)

    #Load up a list of isotopes I want to calculate
    istplist = [
        # 40070, 
        # 60100,
        90180, 
        220450,
        # 521180, 
        # 822000,
        # 471030
        ]

    #Load up the BetaEngine, and pass the isotope list above into the BetaEngine
    BetaSpectraDB = BetaPlusEngine(istplist, xbins = e)

    #Next, I will go ahead and Calculate the BetaSpectrum, and plot out the result, 
    #As well as the total spectrum
    BetaSpectraDB.CalcBetaSpectra(nu_spectrum=True)

    #Iterate over the isotopic list inside the BetaEngine we intitialized, and plot out
    #the Beta Spectrum for each isotope.
    fig = plt.plot()
    for i in istplist:
        betaIstp = BetaSpectraDB.istplist[i]
        plt.errorbar(e, betaIstp.spectrum, betaIstp.uncertainty, label = betaIstp.name)

    plt.xlabel("E (MeV)")
    plt.ylabel("neutrino/MeV")
    plt.legend()
    # plt.show()
    plt.savefig("Individual_bp_neutrino_spectrum.pdf")
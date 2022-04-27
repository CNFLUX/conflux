from conflux.BetaEngine import BetaEngine
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    #I am creating a list of isotopes that I will parse into the
    #Beta Engine. Here, I am using the format ZAI, where
    #the first two numbers will be the atomic number (Z)
    #The second three numbers will be the mass number (A),
    #And the last number is the ____ (I)

    #This list contains Y-96, Te-134, and I-134
    isoList = [390960, 390961, 521340, 531340, 160430, 180500]

    #Here, I am initializing the Beta Engine with my list of isotopes.
    #Note that this is only an example for this specific engine. If you go to
    #the FullReacSpec_SumEngine_Example, you'll see that we initialize the
    #BetaEngine with a specific set of isotopes based off our fission yield engine,
    #whose isotopes we can find in the fissionDB folder.
    Engine = BetaEngine(isoList)
    #Engine.LoadBetaDB()

    #You'll notice, if you've looked in the BetaEngine source file, that I haven't
    #specifically used a function to load the Beta Spectrum database. That's because
    #I call the LoadBetaDB function at the beginning of CalcBetaSpectra.
    #Please also note at this time, you cannot change the binning/threshold of the
    #BetaSpectrum without going directly into the betaEngine source file and manually changing
    #It there, this is so that when we call the summation engine, the list sizes and
    #values line up correctly.
    #Engine.CalcBetaSpectra(targetDB = None, nu_spectrum=True, binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0)
    Engine.CalcBetaSpectra()

    #Now that we've calculated the Beta Spectrum, we can plot our spectrum
    #Note that we've saved our spectra in a dictionary, and to access the spectrum
    #We need to know the specific 'key' for each item. Thankfully, the 'keys' to our
    #dictionary are the isotopes that we used as the input to initialize the betaEngine
    #The spectrum is saved inside the 'spectralist' dictionary, and the uncertainty
    #information is saved inside the 'uncertaintylist'.


    #Here, I'm going to plot out the Beta Spectrum for I-134 and save it as Te-134.png
    fig = plt.figure()
    x = np.linspace(0,20,200)
    plt.errorbar(x, Engine.spectralist[180500], yerr=Engine.uncertaintylist[180500])
    plt.errorbar(x, Engine.spectralist[160430], yerr=Engine.uncertaintylist[160430])

    fig.savefig("Te-134.png")

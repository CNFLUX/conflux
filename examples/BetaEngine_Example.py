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
    isoList = {390960:'Y-96', 521340:'Te-134', 531340:'I-134'}

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
    #You can also change up the range of branches you want to look at. Here, I've
    #set the range from 0MeV to 10MeV
    #Engine.CalcBetaSpectra(self, targetDB = None, nu_spectrum=True,branchErange=[0.0, 20.0])
    Engine.CalcBetaSpectra(nu_spectrum=False, branchErange=[0.0, 10.0])
    print(Engine.istplist)
    #Now that we've calculated the Beta Spectrum, we can plot our spectrum
    #Note that we've saved our spectra in a dictionary, and to access the spectrum
    #We need to know the specific 'key' for each item. Thankfully, the 'keys' to our
    #dictionary are the isotopes that we used as the input to initialize the betaEngine
    #The spectrum is saved inside the 'spectrum' dictionary, and the spectrum uncertainty
    #information is saved inside the 'spectUnc'.


    #Here, I'm going to plot out the Beta Spectrum for I-134 and save it as Te-134.png

    isotope = 531340
    isoname = isoList[isotope]
    fig = plt.figure()
    x = np.linspace(0,20,200)
    #plt.errorbar(x, Engine.spectralist[521340], yerr=Engine.spectUncList[521340])
    # print(Engine.istplist[isotope].spectrum)
    # print(Engine.istplist[isotope].spectUnc)
    plt.plot(x, Engine.istplist[isotope].spectrum)
    plt.fill_between(x, Engine.istplist[isotope].spectrum+Engine.istplist[isotope].spectUnc, Engine.istplist[isotope].spectrum-Engine.istplist[isotope].spectUnc, alpha=.5, linewidth=0)
    fig.savefig(isoname+"s.png")
    fig = plt.figure()
    plt.fill_between(x, Engine.istplist[isotope].spectrum+Engine.istplist[isotope].branchUnc, Engine.istplist[isotope].spectrum-Engine.istplist[isotope].branchUnc, alpha=.5, linewidth=0)
    fig.savefig(isoname+"b.png")
    fig = plt.figure()
    plt.fill_between(x, Engine.istplist[isotope].spectrum+Engine.istplist[isotope].totalUnc, Engine.istplist[isotope].spectrum-Engine.istplist[isotope].totalUnc, alpha=.5, linewidth=0)
    fig.savefig(isoname+"u.png")

    #Finally, some last minute information. If you did happen to print out the istplist inside
    # The BetaEngine, you'll have noticed that there are 4 entries instead of 3. This is
    #Because I-134 has an additional isomer, which is noted in istplist as 531431 instead
    #of 531430
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
if __name__  == "__main__":

    #---------------------------------------------------
    #Full Reactor Spectrum Using U235

    #initialize the Isotopes that you would like to use in your
    #Reaction and load them into the Database
    U235 = FissionIstp(92,235)
    U235.LoadFissionDB()
    #Load up the Correlation data, and calculate the Covariance matrix (both are done by calling LoadCorrelation).
    #Must be done after loading up the Fission Data
    #U235.LoadCorrelation()

    #Initialize the model you'd like to work with, and add the isotope
    #into the model
    FisModel = FissionModel()
    FisModel.AddContribution(isotope=U235, Ei=0, fraction=1)

    #Initialize the type of engine you want to run (The Summation engine in our case)
    #And then add the specific model to that engine
    SummationEngine = SumEngine(xbins=np.arange(0, 15, 0.1))
    SummationEngine.AddModel(FisModel)
    #Before Calculating anything, Normalize the fission products 
    SummationEngine.NormalizeFP()

    #Load in the Beta Shape data and calculate the total beta shape of our reactor.
    #Here I'm inputting a list of isotopes in the FPZAI format from the summation engine
    betaSpectraDB = BetaEngine(SummationEngine.FPYlist.keys(), xbins=np.arange(0, 15, 0.1))
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0,15])

    #Calculate the total reactor spectrum and associated uncertainties
    #from the loaded shape data(betaSpectraDB) And the fission Yield data (SummationEngine)
    SummationEngine.CalcReactorSpectrum(betaSpectraDB, branchErange=[0,15.0])

    #Here I write out the energy bins, spectrum, and uncertainty to a file
    file = open("U235_synth_data.csv", "w")
    file.write("E,Ne,dNe\n")
    for i in range(len(SummationEngine.spectrum)):
        file.write(str(SummationEngine.xbins[i]) + "," + str(SummationEngine.spectrum[i]) + "," + 
                   str(SummationEngine.uncertainty[i]) + "\n")
    file.close()
    
    #This is just a proxy placeholder for the legend << Not necessary for code running
    proxy_1 = Line2D([0], [0], color = 'black', linestyle="-")
    proxy_2 = Line2D([0], [0], color = "blue", linestyle= '--')
    #Draw the resulting spectra
    fig = plt.figure()
    x = SummationEngine.xbins
    for FPZAI in SummationEngine.betaSpectraList:
        #This is to draw out the individual branches <<
        plt.plot(x, SummationEngine.betaSpectraList[FPZAI] * SummationEngine.FPYlist[FPZAI].y)
    #This is the total spectrum
    plt.plot(x, SummationEngine.spectrum, 'b--')
    plt.yscale('log')
    plt.xlabel("E (in MeV)", fontsize= 18)
    plt.ylabel("neutrinos/MeV/Fission", fontsize = 18)
    plt.legend()
    plt.legend(handles=[proxy_1, proxy_2], labels=['Beta Branches', 'Total Spectrum'])
    plt.ylim(bottom= float(10e-13), top = 5)
#    plt.title("U235 Neutrino Spectrum", fontsize=20)
    fig.savefig("U235_Neutrino_Spectrum.png")

    plt.clf()

    #Next, we will pull out all uncertainties and see how they affect our spectrum
    spectrum = SummationEngine.spectrum
    modelUnc = SummationEngine.modelUnc
    totalUnc = SummationEngine.uncertainty
    yieldUnc = SummationEngine.yieldUnc
    plt.plot(x, spectrum, linestyle='-')
    plt.fill_between(x,spectrum + totalUnc, spectrum - totalUnc, label = "totalUnc", alpha = 0.25)
    plt.fill_between(x, spectrum + modelUnc, spectrum - modelUnc, label = "modelUnc", alpha = 0.25)
    plt.fill_between(x, spectrum + yieldUnc, spectrum - yieldUnc, label = "yieldUnc", alpha = 0.25)
    plt.legend()
    plt.savefig("U235_spectrum_with_uncertainties.png")
    
    #Finally, print out all the beta branches that didn't go into creating this spectrum.
    #Each branch is represented in their FPZAI form.
    for i in SummationEngine.missingBranch:
        print(i)

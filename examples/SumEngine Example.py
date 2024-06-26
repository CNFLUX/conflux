from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
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
    model = FissionModel()
    model.AddContribution(isotope=U235, Ei=0, fraction=1)

    #Initialize the type of engine you want to run (The Summation engine in our case)
    #And then add the specific model to that engine
    result = SumEngine(xbins=np.arange(0, 15., 0.1))
    result.AddModel(model)

    #Load in the Beta Shape data and calculate the total beta shape of our reactor
    betaSpectraDB = BetaEngine(result.FPYlist.keys(), xbins=np.arange(0, 15, 0.1))
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 15.0])

    #Calculate the total reactor spectrum from the loaded shape data(betaSpectraDB)
    #And the fission Yield data (result)
    result.CalcReactorSpectrum(betaSpectraDB)

    #Draw the resulting spectra
    fig = plt.figure()
    x = result.xbins
    for FPZAI in result.betaSpectraList:
        #This is to draw out the individual branches
        plt.plot(x, result.betaSpectraList[FPZAI] * result.FPYlist[FPZAI].y)
    #This is the total spectrum
    plt.plot(x, result.spectrum, 'b--')
    plt.yscale('log')
    plt.xlabel("E (in MeV)", fontsize= 18)
    plt.ylabel("neutrinos/MeV/Fission", fontsize = 18)
    plt.title("U235 Neutrino Spectrum", fontsize=20)
    fig.savefig("U235 - Neutrino Spectrum")


    #Finally, print out all the beta branches that didn't go into creating this spectrum.
    #Each branch is represented in their FPZAI form.
    for i in result.missingBranch:
        print(i)

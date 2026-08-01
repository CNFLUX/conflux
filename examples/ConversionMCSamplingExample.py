from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Here, I want to show off the MC sampling aspect of the Conversion example.
    #The more MC samples we have, the better our uncertainty calculation becomes.
    #Most of the code here is the same as the ConversionEngine_Example.py file,
    #So I will not be going in depth into the parts that are the same. For more questions,
    #Please look at that file.

    #Load up the energy range for the calculation (0 to 10 MeV with 100keV steps)
    e = np.linspace(0., 10., 100)


    #Load the beta information from the ConversionDB directory
    #(change this to the local directory for U_235_e_2014.csv)
    beta235 = BetaData("../data/conversionDB/U_235_e_2014.csv")

    #Initialize a U235 Fission Isotope with incident neutron energy of 0 (Thermal)
    U235 = FissionIstp(92,235, Ei=0)
    U235.LoadFissionDB()

    #Load up the conversion model, add the beta data as well as the fission isotpe to
    #The model. the fractional contribution of U235 in this case is 1.
    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)


    #These are some arrays to store the beta and neutrino spectrum information
    betaspec = np.zeros(100)
    nuspec = np.zeros(100)

    #Run conversion mode, incremement the MC sampling for errors by 5 every iteration.
    for i in range(1,10):
        convertmodel.VBfitbeta("U235")
        #here, I increment the number of MC samples I want to take
        sample_size = i * 5
        #Calculate the covariance matrices
        #And here I implement it in the calculation. Note the samples parameter, which will
        #adjust how many MC samples you want to take.
        covmat = convertmodel.vblist["U235"].Covariance(beta235, e, nu_spectrum=False, samples=sample_size)
        covmat_nu = convertmodel.vblist["U235"].Covariance(beta235, e, nu_spectrum=True, samples=sample_size)

        #Here, I calculate the relative beta error by square rooting the leading diagonal terms of my covariance matrix
        #And dividing by the total spectrum
        spectrum = convertmodel.vblist["U235"].SumBranches(e, nu_spectrum=False)
        relativeErr = np.zeros(len(e))
        for j in range(len(spectrum)):
            if spectrum[j] > 0:
                relativeErr[j] = np.sqrt(covmat[j][j]) / spectrum[j]
                betaspec[j] = betaspec[j] + relativeErr[j]

        #Here, I calculate the relative neutrino error, in much the same way I calculated the beta errors.
        spectrumNu = convertmodel.vblist["U235"].SumBranches(e, nu_spectrum=True)
        relativeErrNu = np.zeros(len(e))
        for j in range(len(spectrum)):
            if spectrumNu[j] > 0:
                relativeErrNu[j] = np.sqrt(covmat_nu[j][j]) / spectrumNu[j]
                nuspec[j] = nuspec[j] + relativeErrNu[j]


        #here I plot the individual MC sampling errors for each sample step.
        fig = plt.figure()
        plt.title("MC sample size: " + str(sample_size))
        plt.ylim([-.3, .3])
        plt.fill_between(e, relativeErr, -relativeErr, label='beta', alpha=0.4)
        plt.fill_between(e, relativeErrNu, -relativeErrNu, label='neutrino', alpha=0.4)
        plt.legend()
        name = "U235_MC sample:" + str(i * 5) + "errors.png"
        fig.savefig(name)

        plt.close()
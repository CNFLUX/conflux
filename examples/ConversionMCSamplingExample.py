from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Setup of the conversion mode information (change this to the local directory for U_235_e_2014.csv)
    beta235 = BetaData("../data/conversionDB/U_235_e_2014.csv")

    U235 = FissionIstp(92,235)
    U235.LoadFissionDB()

    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)

    betaspec = np.zeros(200)
    nuspec = np.zeros(200)
    xval = np.linspace(0., 10., 200)

    #Run conversion mode, incremement the MC sampling for errors by 5 every iteration.
    for i in range(1,10):
        xval = np.linspace(0., 10., 200)
        convertmodel.VBfit()
        sample_size = i * 5
        #Calculate the covariance matrices
        covmat = convertmodel.vblist["U235"].Covariance(beta235, xval, nu_spectrum=False, samples=sample_size)
        covmat_nu = convertmodel.vblist["U235"].Covariance(beta235, xval, nu_spectrum=True, samples=sample_size)

        #Calculate the relative beta error
        spectrum = convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum=False)
        relativeErr = np.zeros(len(xval))
        for j in range(len(spectrum)):
            if spectrum[j] > 0:
                relativeErr[j] = np.sqrt(covmat[j][j]) / spectrum[j]
                betaspec[j] = betaspec[j] + relativeErr[j]

        #Calculate the relative neutrino error
        spectrumNu = convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum=True)
        relativeErrNu = np.zeros(len(xval))
        for j in range(len(spectrum)):
            if spectrumNu[j] > 0:
                relativeErrNu[j] = np.sqrt(covmat_nu[j][j]) / spectrumNu[j]
                nuspec[j] = nuspec[j] + relativeErrNu[j]


        #Plot the individual MC sampling errors
        fig = plt.figure()
        plt.title("MC sample size: " + str(sample_size))
        plt.ylim([-.3, .3])
        plt.fill_between(xval, relativeErr, -relativeErr, label='beta', alpha=0.4)
        plt.fill_between(xval, relativeErrNu, -relativeErrNu, label='neutrino', alpha=0.4)
        plt.legend()
        name = "U235_MC sample:" + str(i * 5) + "errors.png"
        fig.savefig(name)

        plt.close()


    negnuspec = []
    negbspec = []
    #Average your MC sampling error to plot
    for i in range(len(nuspec)):
        nuspec[i] = nuspec[i]/5.0
        negnuspec.append(-1. * nuspec[i])
    for j in range(len(betaspec)):
        betaspec[j] = betaspec[j]/5.0
        negbspec.append(-1. * betaspec[j])


    #Plot your averaged MC sampling error
    fig = plt.figure()
    plt.title("Averaged sample")
    plt.ylim([-.3, .3])
    plt.fill_between(xval, nuspec, negnuspec, label='neutrino', alpha=0.4)
    plt.fill_between(xval, betaspec, negbspec, label='beta', alpha=0.4)
    plt.legend()
    name= "averaged sample"
    fig.savefig(name)
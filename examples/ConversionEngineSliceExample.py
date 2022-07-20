from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Load up data to run the conversion engine (change this to the local directory for U_235_e_2014.csv)
    beta235 = BetaData("../data/conversionDB/U_235_e_2014.csv")

    U235 = FissionIstp(92,235)
    U235.LoadFissionDB()

    #Run the conversion engine for U235
    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)

    #Run the virtual beta branch fitting for incrementing branch slices.
    #Branch slices are increased by 0.2, up to 4.0
    for j in range(1,20):
        fig = plt.figure()
        xval = np.linspace(0., 10., 200)
        plt.yscale('log')
        plt.xlabel("E (MeV)")
        plt.ylabel("(electron/neutrino)/fission/MeV")
        k = j * 0.2
        convertmodel.VBfit(k)
        for i in range(0, 15):
            if not sum(convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = False))>0:
                continue
            plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = False), fmt='--')
        plt.errorbar(convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr, label='beta data')
        plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = False), label='beta')
        plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True), label='neutrino')
        plt.legend()
        name = "U235_slice size:" + str(k) + ".png"
        plt.title(name)
        fig.savefig(name)

        covmat = convertmodel.vblist["U235"].Covariance(beta235, xval, nu_spectrum=False, samples=50)
        covmat_nu = convertmodel.vblist["U235"].Covariance(beta235, xval, nu_spectrum=True, samples=50)
        spectrum = convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum=False)
        relativeErr = np.zeros(len(xval))
        for i in range(len(spectrum)):
            if spectrum[i] > 0:
                relativeErr[i] = np.sqrt(covmat[i][i]) / spectrum[i]

        spectrumNu = convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum=True)
        relativeErrNu = np.zeros(len(xval))
        for i in range(len(spectrum)):
            if spectrumNu[i] > 0:
                relativeErrNu[i] = np.sqrt(covmat_nu[i][i]) / spectrumNu[i]

        fig = plt.figure()
        plt.ylim([-.3, .3])
        plt.fill_between(xval, relativeErr, -relativeErr, label='beta', alpha=0.4)
        plt.fill_between(xval, relativeErrNu, -relativeErrNu, label='neutrino', alpha=0.4)
        plt.legend()
        name = "U235_slice size:" + str(k) + "errors.png"
        fig.savefig(name)

        plt.close()
from conflux.FPYEngine import FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData, VirtualBranch
import matplotlib.pyplot as plt
import numpy as np



if __name__ == "__main__":

    #So, we are now in conversion mode. First thing we will have to do is import our
    #Beta data. Thankfully, we've written a method to make this easy, simply
    #specify the path where the data is located in the function.
    #betaU235 = BetaData("path/to/file")
    #(change this to the local directory for U_235_e_2014.csv)
    betaU235 = BetaData("../data/conversionDB/U_235_e_2014.csv")

    #Next, just to make sure the data looks good, we're going to print out the
    #error and the yield.

    # for i in range(len(betaU235.y)):
    #     print("Yield:",betaU235.y[i],"Error:",  betaU235.yerr[i])

    #Now, we'll load up the Fission isotope data for U-235.

    U235 = FissionIstp(92,235)
    U235.LoadFissionDB()

    print("Check One")

    #Now I'm going to initialize the conversion Engine, and add the U235 Fission
    #and beta information into it. Note that the first variable needs to be the
    #beta data, and the second variable needs to be the fission information

    convertmodel = ConversionEngine()
    #convertModel.AddBetaData(BetaData, FissionIstp, str name, int frac)
    convertmodel.AddBetaData(betaU235, U235, "U235", 1.0)

    #Finally, we fit virtual beta branches to the beta spectra data with a user
    #Defined slice size (Default slice size is 0.5). Here I've chosen a slice size of
    #1.0
    xval = np.arange(0.,10.0, 0.05)


    convertmodel.VBfitbeta("U235")

    print("Check Three")
    # Everything above here is the conversion  mode calculation. Once I've done the
    # Calculation, I need to extract the data and plot it/visualize/represent it.
    # All plotting and data-visualization will be below.

    fig = plt.figure()



    #Below I'm going to plot out the total neutrino and beta spectrum, as well as
    #print out all the virtual branches as well. I've also gone ahead and plotted
    #Out the beta data as well.
    plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")
    #Plot out the model branches


    #Something wrong with the plotting on this one


    for i in range(0, 20):
        if not sum(convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = False))>0:
            continue
        plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = False), fmt='--')
    #Plot out the raw beta data
    plt.errorbar(convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr, label='beta data')
    #Plot out the calculated beta spectrum
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = False), label='beta')
    #Plot out the calculated neutrino spectrum
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True), label='neutrino')
    plt.legend()
    fig.savefig("U235_conversion.png")

    #This next set of plots are going to be the covariance matrices for both
    #The neutrino and beta spectrum for U235. I've initialized them here, I
    #will plot them later in the code.

    #beta covariance matrix
    covmat_beta = convertmodel.vblist["U235"].Covariance(betaU235, xval, nu_spectrum=False, samples=50)
    #neutrino covariance matrix
    covmat_nu = convertmodel.vblist["U235"].Covariance(betaU235, xval, nu_spectrum=True, samples=50)


    spectrum = convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum=False)
    relativeErr = np.zeros(len(xval))
    for i in range(len(xval)):
        if spectrum[i] >0:
            relativeErr[i] = np.sqrt(covmat_beta[i][i]/spectrum[i])

    spectrumNu = convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum=True)
    relativeErrNu = np.zeros(len(xval))
    for i in range(len(xval)):
        if spectrumNu[i] > 0:
            relativeErrNu[i] = np.sqrt(covmat_nu[i][i] / spectrumNu[i])

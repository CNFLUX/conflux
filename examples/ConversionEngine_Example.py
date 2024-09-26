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

    #I am also going to intialize my energy range, in this case 0 MeV to 10 MeV with 
    #100 keV bins

    e = np.arange(0, 10., 0.1)

    #Next, just to make sure the data looks good, we're going to print out the
    #error and the yield.
    # for i in range(len(betaU235.y)):
    #     print("Yield:",betaU235.y[i],"Error:",  betaU235.yerr[i])

    #Now, we'll load up the Fission isotope data for U-235, with an incident neutron energy
    #of 0 (thermal), and load the correlation information as well.
    U235 = FissionIstp(92,235, Ei=0)
    U235.LoadFissionDB()
    U235.LoadCorrelation()

    #Next, I need to calculate the spectral shape for all the products of U235, so I will load
    #up a beta engine and calculate that before moving onto the conversion engine step

    #Now I'm going to initialize the conversion Engine, and add the U235 Fission
    #and beta information into it. Note that the first variable needs to be the
    #beta data, and the second variable needs to be the fission information

    convertmodel = ConversionEngine()
    #convertModel.AddBetaData(BetaData, FissionIstp, str name, int frac)
    convertmodel.AddBetaData(betaU235, U235, "U235", 1.0)

    #Finally, we fit virtual beta branches to the inputted fission isotope with
    #a slice size of 0.5 (500keV, which is the default). Note also, that
    #My parameter is a string of the inputted fission isotope, not the isotope
    #Itself. Be careful when adding this parameter.
    convertmodel.VBfitbeta("U235")

    #Now, since I've only added one fission isotope to my conversion model, I do not
    #Need to call the below function. I can infact calculate the spectrum as I have in
    #The plotting section and be done. However, in calculations that require multiple
    #Fission isotopes to be added, it is better to call SummedSpectra to calculate the total
    #spectrum, uncertainty, and covariance matrices for the entire model, not just an
    #individual Fission isotope
    conSpec, conUnc, conCov = convertmodel.SummedSpectrum(e, nu_spectrum=True, cov_samp=5)


    #Now, I will start plotting.
    fig = plt.figure()

    #Below I'm going to plot out the total neutrino and beta spectrum, as well as
    #print out all the virtual branches as well.
    plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")


    #This for loop is to plot out the individual beta branches of this reaction.
    for i in range(0, 20):
        #This checks to see if the sum of all the virtual branches is positive. If it is not, then skip
        if not sum(convertmodel.vblist["U235"].SumBranches(e, thresh =i*0.5, nu_spectrum = False))>0:
            continue
        #Otherwise, plot out this branch by summing all the virtual beta branches whose energy range is above
        #the threshold.
        plt.errorbar(e, convertmodel.vblist["U235"].SumBranches(e, thresh =i*0.5, nu_spectrum = False),fmt='--')
    
    #Here, I plot out the raw beta data
    plt.errorbar(convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr, label='beta data')
    
    #Plot out the calculated beta spectrum
    plt.errorbar(e, convertmodel.vblist["U235"].SumBranches(e, nu_spectrum = False), label='beta')
    #Plot out the calculated neutrino spectrum
    plt.errorbar(e, convertmodel.vblist["U235"].SumBranches(e, nu_spectrum = True), label='neutrino')
    plt.legend()
    fig.savefig("U235_conversion.png")

    #Note, I am going to plot out the converted spectrum that we got from SummedSpectra, and will show that it's the
    #Exact same as if I called SumBranches from the last plt.errorbar call above.

    plt.clf()
    plt.errorbar(e, conSpec, yerr = conUnc, label = "SummedSpec", fmt = "--")
    plt.errorbar(e, convertmodel.vblist["U235"].SumBranches(e, nu_spectrum = True), label='SumBranches')
    plt.legend()
    fig.savefig("comparison.png")


    #This calculates the errors of our conversion calcualtion, with 50 samples.
    #I generate a covariance matrix for both the betas and the neutrinos
    covmat = convertmodel.vblist["U235"].Covariance(betaU235, e, nu_spectrum=False, samples=50)
    covmat_nu = convertmodel.vblist["U235"].Covariance(betaU235, e, nu_spectrum=True, samples=50)

    #Here, I calculate the beta spectrum, and then calculate the relative errors of it from the
    #Covariance matrices I generated above.
    spectrumBeta = convertmodel.vblist["U235"].SumBranches(e, nu_spectrum=False)
    relativeErr = np.zeros(len(e))
    for i in range(len(spectrumBeta)):
        if spectrumBeta[i] > 0:
            relativeErr[i] = np.sqrt(covmat[i][i]) / spectrumBeta[i]


    #Here, I calculate the neutrino spectrum, and then calculate the relative errors of it from the
    #Covariance matrices I generated above.
    spectrumNu = convertmodel.vblist["U235"].SumBranches(e, nu_spectrum=True)
    relativeErrNu = np.zeros(len(e))
    for i in range(len(spectrumNu)):
        if spectrumNu[i] > 0:
            relativeErrNu[i] = np.sqrt(covmat_nu[i][i]) / spectrumNu[i]


    #Finally, the rest of this is plotting the relative errors for both the betas and the neutrinos for our given
    #Energy range
    fig = plt.figure()
    plt.ylim([-.3, .3])
    plt.fill_between(e, relativeErr, -relativeErr, label='beta', alpha=0.4)
    plt.fill_between(e, relativeErrNu, -relativeErrNu, label='neutrino', alpha=0.4)
    plt.legend()
    name = "U235_errors.png"
    fig.savefig(name)

    plt.close()
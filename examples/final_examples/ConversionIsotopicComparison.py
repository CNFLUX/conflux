from conflux.config import CONFLUX_DB
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
    betaU235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")
    betaPu239 = BetaData(CONFLUX_DB+"/conversionDB/Pu_239_e_2014.csv")
    betaPu241 = BetaData(CONFLUX_DB+"/conversionDB/Pu_241_e_2014.csv")

    #I am also going to intialize my energy range, in this case 0 MeV to 10 MeV with 
    #100 keV bins

    e = np.arange(0, 10., 0.1)

    #Next, just to make sure the data looks good, we're going to print out the
    #error and the yield.
    # for i in range(len(betaU235.y)):
    #     print("Yield:",betaU235.y[i],"Error:",  betaU235.yerr[i])

    #Now, we'll load up the Fission isotope data for U-235, with an incident neutron energy
    #of 0 (thermal), and load the correlation information as well.
    U235 = FissionIstp(92, 235, Ei=0)
    U235.LoadFissionDB()
    Pu239 = FissionIstp(94, 239, Ei=0)
    Pu239.LoadFissionDB()
    Pu241 = FissionIstp(94, 241, Ei=0)
    Pu241.LoadFissionDB()

    #Next, I need to calculate the spectral shape for all the products of U235, so I will load
    #up a beta engine and calculate that before moving onto the conversion engine step

    #Now I'm going to initialize the conversion Engine, and add the U235 Fission
    #and beta information into it. Note that the first variable needs to be the
    #beta data, and the second variable needs to be the fission information

    convertmodel = ConversionEngine()
    #convertModel.AddBetaData(BetaData, FissionIstp, str name, int frac)
    convertmodel.AddBetaData(betaU235, U235, "U235", 1.0)
    convertmodel.AddBetaData(betaPu239, Pu239, "Pu239", 1.0)
    convertmodel.AddBetaData(betaPu241, Pu241, "Pu241", 1.0)


    #Finally, we fit virtual beta branches to the inputted fission isotope with
    #a slice size of 0.5 (500keV, which is the default). Note also, that
    #My parameter is a string of the inputted fission isotope, not the isotope
    #Itself. Be careful when adding this parameter.
    convertmodel.VBfitbeta("U235", slicesize=0.5)
    convertmodel.VBfitbeta("Pu239", slicesize=0.5)
    convertmodel.VBfitbeta("Pu241", slicesize=0.5)
    

    #Now, since I've only added one fission isotope to my conversion model, I do not
    #Need to call the below function. I can infact calculate the spectrum as I have in
    #The plotting section and be done. However, in calculations that require multiple
    #Fission isotopes to be added, it is better to call SummedSpectra to calculate the total
    #spectrum, uncertainty, and covariance matrices for the entire model, not just an
    #individual Fission isotope
    conSpec, conUnc, conCov = convertmodel.SummedSpectrum(e, nu_spectrum=False, cov_samp=1)

    fig = plt.figure()    

    plt.errorbar(betaU235.xbins, betaU235.spectrum, yerr=betaU235.uncertainty, label='U235')
    plt.errorbar(betaPu239.xbins, betaPu239.spectrum, yerr=betaPu239.uncertainty, label='Pu239')
    plt.errorbar(betaPu241.xbins, betaPu241.spectrum, yerr=betaPu241.uncertainty, label='Pu239')
    plt.legend()
    plt.xlabel("E (MeV)")
    plt.ylabel("beta/fission/MeV")
    plt.show()
    fig.savefig("betaspectra.pdf")

    fig = plt.figure()    

    plt.errorbar(e, convertmodel.vblist["U235"].spectrum, label='U235')
    plt.errorbar(e, convertmodel.vblist["Pu239"].spectrum, label='Pu239')
    plt.errorbar(e, convertmodel.vblist["Pu241"].spectrum, label='Pu241')
    plt.legend()
    # plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.ylabel("neutrino/fission/MeV")
    fig.savefig("isotopic_comparison.pdf")

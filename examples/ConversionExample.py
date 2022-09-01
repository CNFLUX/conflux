# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

# test
if __name__ == "__main__":
    # Begin the calculation by sourcing the default beta data
    beta235 = BetaData("./data/conversionDB/U_235_e_2014.csv")
    beta239 = BetaData("./data/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData("./data/conversionDB/Pu_241_e_2014.csv")
    
    # Define isotopic fission yield DB to calculate average atom numbers of
    # virtual branches
    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)
    U235.LoadFissionDB()
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()

    # Declare the conversion engine by adding beta data with corresponding FPY
    # database
    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta239, Pu239, "Pu239", 1.0)
    convertmodel.VBfit(0.5)

    # Draw plots to test output.
    print("Drawing spectrum...")
    fig = plt.figure()
    xval = np.linspace(0., 10., 200)

    plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")
    for i in range(0, 40):
        if (sum(convertmodel.vblist["Pu239"].SumBranches(xval,
                thresh =i*0.25, nu_spectrum = False))<=0):
            continue
        plt.errorbar(xval, convertmodel.vblist["Pu239"].SumBranches(xval,
            thresh =i*0.25, nu_spectrum = False), fmt='--')
    plt.errorbar(convertmodel.betadata["Pu239"].x,
        convertmodel.betadata["Pu239"].y, convertmodel.betadata["Pu239"].yerr,
        label='beta data')
    plt.errorbar(xval, convertmodel.vblist["Pu239"].SumBranches(xval,
        nu_spectrum = False), label='beta')
    plt.errorbar(xval, convertmodel.vblist["Pu239"].SumBranches(xval,
        nu_spectrum = True), label='neutrino')
    plt.legend()
    fig.savefig("Pu239_conversion.png")

    # ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    covmat=convertmodel.vblist["Pu239"].Covariance(beta239,
        xval, nu_spectrum = False, samples=50)
    covmat_nu=convertmodel.vblist["Pu239"].Covariance(beta239,
        xval, nu_spectrum = True, samples=50)
    #print(covmat)
    
    spectrum = convertmodel.vblist["Pu239"].SumBranches(xval,
        nu_spectrum = False)
    relativeErr = np.zeros(len(xval))
    for i in range(len(spectrum)):
        if spectrum[i] > 0:
            relativeErr[i] = np.sqrt(covmat[i][i])/spectrum[i]
    
    spectrumNu = convertmodel.vblist["Pu239"].SumBranches(xval,
    nu_spectrum = True)
    relativeErrNu = np.zeros(len(xval))
    for i in range(len(spectrum)):
        if spectrumNu[i] > 0:
            relativeErrNu[i] = np.sqrt(covmat_nu[i][i])/spectrumNu[i]
    
    fig = plt.figure()
    plt.ylim([-.3, .3])
    plt.fill_between(xval, relativeErr, -relativeErr, label='beta',
    alpha = 0.4)
    plt.fill_between(xval, relativeErrNu, -relativeErrNu, label='neutrino', 
    alpha = 0.4)
    plt.legend()
    fig.savefig("Pu239_errs.png")
        
    fig = plt.figure()
    im = plt.imshow(covmat)
    plt.colorbar(im)
    fig.savefig("Pu239_cov.png")
    
    fig = plt.figure()
    im = plt.imshow(covmat_nu)
    plt.colorbar(im)
    fig.savefig("Pu239_cov_nu.png")

# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

# test
if __name__ == "__main__":
    beta235 = BetaData("./data/conversionDB/U_235_e_2014.csv")
    beta239 = BetaData("./data/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData("./data/conversionDB/Pu_241_e_2014.csv")
    #print(beta235.y)
    #print(beta235.yerr)

    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)
    U235.LoadFissionDB()
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()

    # vbtest = VirtualBranch(U235)
    # vbtest.CalcZAavg(6,7)
    # print(vbtest.Aavg, vbtest.Zavg)

    
    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta239, Pu239, "Pu239", 1.0)
    convertmodel.VBfit(0.25)

    # Draw plots to test output.
    print("Drawing spectrum...")
    fig = plt.figure()
    xval = np.linspace(0., 10., 200)
    #plt.errorbar(self.betadata[name].x, self.vblist[name].SumBranches(self.betadata[name].x, True))
    # for i in range(0, 20):
    #     print(self.vblist[name].SumBranches(xval, thresh =i*0.5, nu_spectrum = False))
    #     plt.errorbar(xval, self.vblist[name].SumBranches(xval, thresh =i*0.5, nu_spectrum = False), fmt='--')
    plt.yscale('log')
    plt.ylim([1e-5, 2])
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")
    for i in range(0, 40):
        if not sum(convertmodel.vblist["Pu239"].SumBranches(xval, thresh =i*0.25, nu_spectrum = False))>0:
            continue
        plt.errorbar(xval, convertmodel.vblist["Pu239"].SumBranches(xval, thresh =i*0.25, nu_spectrum = False), fmt='--')
    plt.errorbar(convertmodel.betadata["Pu239"].x, convertmodel.betadata["Pu239"].y, convertmodel.betadata["Pu239"].yerr, label='beta data')
    plt.errorbar(xval, convertmodel.vblist["Pu239"].SumBranches(xval, nu_spectrum = False), label='beta')
    plt.errorbar(xval, convertmodel.vblist["Pu239"].SumBranches(xval, nu_spectrum = True), label='neutrino')
    plt.legend()
    fig.savefig("Pu239_conversion_new.png")

    # ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    covmat=convertmodel.vblist["Pu239"].Covariance(beta239, xval, nu_spectrum = False, samples=50)
    covmat_nu=convertmodel.vblist["Pu239"].Covariance(beta239, xval, nu_spectrum = True, samples=50)
    #print(covmat)
    
    spectrum = convertmodel.vblist["Pu239"].SumBranches(xval, nu_spectrum = False)
    relativeErr = np.zeros(len(xval))
    for i in range(len(spectrum)):
        if spectrum[i] > 0:
            relativeErr[i] = np.sqrt(covmat[i][i])/spectrum[i]
    
    spectrumNu = convertmodel.vblist["Pu239"].SumBranches(xval, nu_spectrum = True)
    relativeErrNu = np.zeros(len(xval))
    for i in range(len(spectrum)):
        if spectrumNu[i] > 0:
            relativeErrNu[i] = np.sqrt(covmat_nu[i][i])/spectrumNu[i]
    
    fig = plt.figure()
    plt.ylim([-.3, .3])
    plt.fill_between(xval, relativeErr, -relativeErr, label='beta', alpha = 0.4)
    plt.fill_between(xval, relativeErrNu, -relativeErrNu, label='neutrino', alpha = 0.4)
    plt.legend()
    fig.savefig("Pu239_errs_new.png")
        
    fig = plt.figure()
    im = plt.imshow(covmat)
    plt.colorbar(im)
    fig.savefig("Pu239_cov_new.png")
    
    fig = plt.figure()
    im = plt.imshow(covmat_nu)
    plt.colorbar(im)
    fig.savefig("Pu239_cov_nu_new.png")

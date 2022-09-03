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
    
    # Loading default fission product DB
    U235.LoadFissionDB()
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()

    # Define the size of energy slice
    branch_slice = 0.25
    # Declare the conversion engine by adding beta data with corresponding FPY
    # database
    convertmodel = ConversionEngine()
    # Add beta spectra and fission products to the conversion engine
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)
    # Do virtual branch fitting with the defined virtual branch energy range
    convertmodel.VBfit(branch_slice)

    # Draw plots to test output.
    print("Drawing spectrum...")
    fig = plt.figure()
    # Define x values of data points to be drawn in the figure
    xval = np.arange(2., 10., 0.25)
    
    # Setting figure parameters
    plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")
    
    # # Draw the spectra of all vertual branches
    # for i in range(0, 40):
    #     if (sum(convertmodel.vblist["Pu239"].SumBranches(xval,
    #             thresh =i*branch_slice, nu_spectrum = False))<=0):
    #         continue
    #     plt.errorbar(xval, convertmodel.vblist["Pu239"].SumBranches(xval,
    #         thresh =i*branch_slice, nu_spectrum = False), fmt='--')
    #
    # # Draw the original beta spectrum
    # plt.errorbar(convertmodel.betadata["Pu239"].x,
    #     convertmodel.betadata["Pu239"].y, convertmodel.betadata["Pu239"].yerr,
    #     label='beta data')
    # # Draw the summed best fit virtual beta spectrum of this calculation
    # spectrum = convertmodel.vblist["Pu239"].SumBranches(xval,
    #     nu_spectrum = False)
    # plt.errorbar(xval, spectrum, label='best fit beta')
    # # Draw the summed neutrino spectrum from converted best fit beta spectra
    # spectrumNu = convertmodel.vblist["Pu239"].SumBranches(xval,
    # nu_spectrum = True)
    # plt.errorbar(xval, spectrumNu, label='neutrino')
    # plt.legend()
    # fig.savefig("Pu239_conversion_new.png")
    #
    # # Calculate covariance matrix of summed best fit beta spectrum
    # covmat=convertmodel.vblist["Pu239"].Covariance(beta239,
    #     xval, nu_spectrum = False, samples=50)
    # # Calculate covariance matrix of summed converted neutrino spectrum
    # covmat_nu=convertmodel.vblist["Pu239"].Covariance(beta239,
    #     xval, nu_spectrum = True, samples=50)
    #
    # relativeErr = np.zeros(len(xval))
    # for i in range(len(spectrum)):
    #     if spectrum[i] > 0:
    #         relativeErr[i] = np.sqrt(covmat[i][i])/spectrum[i]
    #
    # relativeErrNu = np.zeros(len(xval))
    # for i in range(len(spectrum)):
    #     if spectrumNu[i] > 0:
    #         relativeErrNu[i] = np.sqrt(covmat_nu[i][i])/spectrumNu[i]
    #
    # fig = plt.figure()
    # plt.ylim([-.3, .3])
    # plt.fill_between(xval, relativeErr, -relativeErr, label='beta',
    # alpha = 0.4)
    # plt.fill_between(xval, relativeErrNu, -relativeErrNu, label='neutrino',
    # alpha = 0.4)
    # plt.legend()
    # fig.savefig("Pu239_errs_new.png")
    #
    # fig = plt.figure()
    # im = plt.imshow(covmat)
    # plt.colorbar(im)
    # fig.savefig("Pu239_cov_new.png")
    #
    # fig = plt.figure()
    # im = plt.imshow(covmat_nu)
    # plt.colorbar(im)
    # fig.savefig("Pu239_cov_nu_new.png")
    
    final_spect, final_unc, final_cov = convertmodel.SummedSpectrum(xval)
    print(final_spect)
    print(final_unc)
    print(final_cov)

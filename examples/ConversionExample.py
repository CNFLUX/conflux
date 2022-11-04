# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

def HuberZavg(x, c0, c1, c2):
    return c0+c1*x+c2*x**2

def Rebin(inputx, inputy, outputx):
    # Rebinning the histogram
    new_xbins = outputx
    new_bin_width = outputx[1] - outputx[0]
    new_xbins_low= new_xbins - new_bin_width/2
    new_ybins = np.zeros(len(new_xbins))
    
    # get information from the old bins
    old_xbins = inputx
    old_ybins = inputy
    old_bin_width = old_xbins[1] - old_xbins[0]
    old_bins_low = inputx - old_bin_width/2

    for j in range(len(new_xbins_low)):
        norm = 0
        for i in range(len(old_bins_low)):
            # process if only old bin is in the new bin's range
            weight = old_bin_width
            if old_bins_low[i]>=new_xbins_low[j]:
                # deal with case if old bin is partially in the new bin
                old_bin = old_bins_low[i]
                new_bin = new_xbins_low[j]
                diff = old_bin - new_bin
                
                if diff < old_bin_width:
                    weight = diff
                    norm += weight
                    # print(1, i, old_bin, new_bin, weight, old_ybins[i])
                    new_ybins[j] += weight * old_ybins[i]
                    
                elif diff >= new_bin_width:
                    weight = old_bin_width + new_bin_width - diff
                    norm += weight
                    # print(2, i, old_bin, new_bin, weight, old_ybins[i])
                    new_ybins[j] += weight * old_ybins[i]
                    break
                    
                else:
                    norm += weight
                    # print(3, i, old_bin, new_bin, weight, old_ybins[i])
                    new_ybins[j] += weight * old_ybins[i]

        if norm > 0:
            new_ybins[j] /= norm
    
    return new_ybins

# test
if __name__ == "__main__":
    # Begin the calculation by sourcing the default beta data
    beta235 = BetaData("./data/conversionDB/U_235_e_2014.csv")
    beta235s = BetaData("./data/conversionDB/Synthetic_235_beta.csv")
    beta239 = BetaData("./data/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData("./data/conversionDB/Pu_241_e_2014.csv")
    
    # Define isotopic fission yield DB to calculate average atom numbers of
    # virtual branches
    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)
    
    # Loading default fission product DB
    U235.LoadFissionDB(defaultDB='JEFF')
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
    xval = np.arange(0,10, 0.01)
    Zlist = dict(zip(xval, HuberZavg(xval, 49, -0.4, -0.084)))
    convertmodel.VBfitbeta("U235", branch_slice, Zlist=Zlist)

    # Draw plots to test output.
    print("Drawing spectrum...")
    # Define x values of data points to be drawn in the figure
    
    
    # Setting figure parameters
    fig = plt.figure()
    # plt.yscale('log')
    # plt.ylim([1e-5, 2])
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")
    
    # Draw the spectra of all vertual branches
    for i in range(0, 40):
        if (sum(convertmodel.vblist["U235"].SumBranches(xval,
                thresh=i*branch_slice, nu_spectrum = False))<=0):
            continue
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
            thresh=1*branch_slice, nu_spectrum = False), fmt='--')
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
            thresh=1*branch_slice, nu_spectrum = True), fmt='--')
    
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
    # final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval)

    final_spect, final_unc, final_cov = convertmodel.SummedSpectrum(xval, nu_spectrum=False, cov_samp=20)
    # final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval, nu_spectrum=True)

    # print(final_spect)
    # print(final_spect)
    
    # print(final_unc)
    # print(final_cov)
    print("beta integral", sum(final_spect))
    # print("neu integral",sum(final_spect1))

    
    newxval = np.arange(2, 8.25, 0.25)
    newyval = Rebin(xval, final_spect, newxval)
    # for i in convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True):
    #     print(i)
    
    # testxval = np.linspace(0,200,201)
    # testyval = 1*testxval+2
    # testxout = np.linspace(0,200,21)
    # testyout = Rebin(testxval, testyval, testxout)
    
    fig = plt.figure()
    # plt.yscale('log')
    plt.errorbar(convertmodel.betadata["U235"].x,
        convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
        label='beta data')
    plt.plot(newxval, newyval)
    # plt.plot(xval, final_spect1)
    plt.show()
    
    betaspect = np.interp(xval, convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y)
    diff = (final_spect-betaspect)/betaspect
    fig = plt.figure()
    plt.ylim([-0.2, 0.2])
    plt.plot(xval, diff)
    plt.show()

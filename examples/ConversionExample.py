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
    beta235 = BetaData("../../data/conversionDB/U_235_e_2014.csv")
    beta235 = BetaData("./U235_synth_data_50keV.csv")
    beta239 = BetaData("../../data/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData("../../data/conversionDB/Pu_241_e_2014.csv")
    
    # Define isotopic fission yield DB to calculate average atom numbers of
    # virtual branches
    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)
    
    # Loading default fission product DB
    U235.LoadFissionDB()
    U235.LoadCorrelation()
    U235.CalcCovariance(Ei=0)
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()

    # Define the size of energy slice
    branch_slice = 0.25
    # Declare the conversion engine by adding beta data with corresponding FPY
    # database
    convertmodel = ConversionEngine()
    # Add beta spectra and fission products to the conversion engine
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)
    # convertmodel.AddBetaData(beta239, Pu239, "Pu239", 1.0)
    # convertmodel.AddBetaData(beta241, Pu241, "Pu241", 1.0)
    # Do virtual branch fitting with the defined virtual branch energy range
    xval = np.arange(0,10, 0.05)
    Zlist = dict(zip(xval, HuberZavg(xval, 49, -0.4, -0.084)))
    #convertmodel.VBfitbeta("U235", branch_slice)
    convertmodel.VBfitbeta("U235", branch_slice)


    # Draw plots to test output.
    print("Drawing spectrum...")
    # Define x values of data points to be drawn in the figure
    
    
    # Setting figure parameters
    fig = plt.figure()
    # plt.yscale('log')
    plt.xlim([2, 10])
    plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")
    
    # # Draw the spectra of all vertual branches
    # for i in range(0, 40):
    #     if (sum(convertmodel.vblist["U235"].SumBranches(xval,
    #             thresh=i*branch_slice, nu_spectrum = False))<=0):
    #         continue
    # plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
    #         thresh=1*branch_slice, nu_spectrum = False), fmt='--')
    # plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
    #         thresh=1*branch_slice, nu_spectrum = True), fmt='--')
    #
    # # Draw the original beta spectrum
    # plt.errorbar(convertmodel.betadata["U235"].x,
    #     convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
    #     label='beta data')
    # # Draw the summed best fit virtual beta spectrum of this calculation
    # spectrum = convertmodel.vblist["U235"].SumBranches(xval,
    #     nu_spectrum = False)
    # plt.errorbar(xval, spectrum, label='best fit beta')
    # # Draw the summed neutrino spectrum from converted best fit beta spectra
    # spectrumNu = convertmodel.vblist["U235"].SumBranches(xval,
    # nu_spectrum = True)
    # plt.errorbar(xval, spectrumNu, label='neutrino')
    # plt.legend()
    # fig.savefig("U235_conversion_new.png")



    # fig = plt.figure()
    # im = plt.imshow(covmat)
    # plt.colorbar(im)
    # fig.savefig("U235_cov_new.png")
    #
    # fig = plt.figure()
    # im = plt.imshow(covmat_nu)
    # plt.colorbar(im)
    # fig.savefig("U235_cov_nu_new.png")
    # final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval)


    #Pull the neutrino and beta spectrum out of the conversion engine.
    final_spect, final_unc, final_cov = convertmodel.SummedSpectrum(xval, nu_spectrum=False, cov_samp=25)
    final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval, nu_spectrum=True, cov_samp=25)

    #Print out both spectra
    print(final_spect)
    print(final_cov1)

    # print(final_spect)
    # print(final_spect)
    # Calculate covariance matrix of summed best fit beta spectrum
    covmat=convertmodel.vblist["U235"].Covariance(beta235,
        xval, nu_spectrum = False, samples=25)
    # Calculate covariance matrix of summed converted neutrino spectrum
    covmat_nu=convertmodel.vblist["U235"].Covariance(beta235,
        xval, nu_spectrum = True, samples=25)
    #
    #Calculate the relative error for both the neutrino and beta spectrum, and then
    #Plot both out.
    #


    relativeErr = np.zeros(len(xval))
    for i in range(len(final_spect)):
        if final_spect[i] > 0:
            relativeErr[i] = np.sqrt(final_cov[i][i])/final_spect[i]
    
    relativeErrNu = np.zeros(len(xval))
    for i in range(len(final_spect1)):
        if final_spect1[i] > 0:
            relativeErrNu[i] = np.sqrt(final_cov1[i][i])/final_spect1[i]
    
    fig = plt.figure()
    plt.ylim([-.3, .3])
    plt.fill_between(xval, relativeErr, -relativeErr, label='beta',
                    alpha = 0.4)
    plt.fill_between(xval, relativeErrNu, -relativeErrNu, label='neutrino',
                    alpha = 0.4)
    plt.legend()
    fig.savefig("U235_errs_new.png")
    # print(final_unc)
    # print(final_cov)
    print("beta integral", sum(final_spect))
    # print("neu integral",sum(final_spect1))

    
    newxval = np.arange(2.125, 8.375, 0.25)
    newyval = Rebin(xval, final_spect, newxval)
    newyval1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(newxval, nu_spectrum=True, cov_samp=20)
    # for i in convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True):
    #     print(i)
    
    # testxval = np.linspace(0,200,201)
    # testyval = 1*testxval+2
    # testxout = np.linspace(0,200,21)
    # testyout = Rebin(testxval, testyval, testxout)
    for a in (newyval1):
        print(a)
    
    fig = plt.figure()
    # plt.yscale('log')
    plt.errorbar(convertmodel.betadata["U235"].x,
        convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
        label='beta data')
    plt.plot(newxval, newyval)
    # plt.plot(xval, final_spect1)
    plt.savefig("CE2.png")
    
    betaspect = np.interp(xval, convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y)
    diff = (final_spect-betaspect)/betaspect
    fig = plt.figure()
    plt.ylim([-0.05, 0.05])
    plt.xlim([2, 9])
    plt.plot(xval, diff)
    plt.savefig("CE1.png")

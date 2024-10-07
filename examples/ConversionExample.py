# local modules
from conflux.BetaEngine import BetaEngine, CONFLUX_DB, BetaBranch
from conflux.FPYEngine import FissionIstp
from conflux.SumEngine import SumEngine
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
    beta235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")
    beta235s = BetaData("./U235_synth_data_1.5_9.6.csv")
    beta239 = BetaData(CONFLUX_DB+"/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData(CONFLUX_DB+"/conversionDB/Pu_241_e_2014.csv")

    # Define isotopic fission yield DB to calculate average atom numbers of
    # virtual branches
    U235 = FissionIstp(92, 235, Ei=0)
    Pu239 = FissionIstp(94, 239, Ei=0)
    Pu241 = FissionIstp(94, 241, Ei=0)

    # Loading default fission product DB
    U235.LoadFissionDB(DB='JEFF')
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()

    # Define the size of energy slice
    branch_slice = 0.25
    # Declare the conversion engine by adding beta data with corresponding FPY
    # database
    convertmodel = ConversionEngine()
    # Add beta spectra and fission products to the conversion engine
    convertmodel.AddBetaData(beta235s, U235, "U235", 1.0)
    # convertmodel.AddBetaData(beta239, Pu239, "Pu239", 1.0)
    # convertmodel.AddBetaData(beta241, Pu241, "Pu241", 1.0)
    # Do virtual branch fitting with the defined virtual branch energy range
    xval = np.arange(0, 10, 0.02)
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
    plt.ylim([1e-5, 10])
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")

    # Draw the spectra of all vertual branches
    vbfitter = convertmodel.vblist["U235"]
    vbspec = np.zeros(len(xval))
    binwidths = xval[1]-xval[0]
    for s in vbfitter.E0:
        vbcont = vbfitter.contribute[s]
        vbranch = BetaBranch(vbfitter.Zlist[s], 
                                     vbfitter.Alist[s],
                                     frac=1-vbfitter.fblist[s], 
                                     I=0, 
                                     Q = vbfitter.E0[s],
                                     E0=vbfitter.E0[s], 
                                     sigma_E0=0, 
                                     sigma_frac=0,
                                     forbiddenness=0, 
                                     bAc=vbfitter.wmlist[s])
        spec = vbranch.BetaSpectrum(xval, nu_spectrum=False)
        spec = spec/(sum(spec)*binwidths)
        vbspec += vbcont * spec
        plt.plot(xval, vbspec, "--")
        
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
            thresh=1*branch_slice, nu_spectrum = False), label = "best-fit beta")
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
            thresh=1*branch_slice, nu_spectrum = True), label = "converted_neutrino")
    fig.savefig("U235_convert_test.png")
    
    fig = plt.figure()
    # Draw the original beta spectrum
    plt.errorbar(convertmodel.betadata["U235"].x,
        convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
        label='beta data')
    
    # Draw the summed best fit virtual beta spectrum of this calculation 
    spectrum = convertmodel.vblist["U235"].SumBranches(xval,
        nu_spectrum = False)
    plt.errorbar(xval, spectrum, label='best fit beta')
    # Draw the summed neutrino spectrum from converted best fit beta spectra
    spectrumNu = convertmodel.vblist["U235"].SumBranches(xval,
    nu_spectrum = True)
    plt.errorbar(xval, spectrumNu, label='neutrino')
    plt.legend()
    fig.savefig("U235_conversion_new.png")

    final_spect, final_unc, final_cov = convertmodel.SummedSpectrum(xval, 
                                        nu_spectrum=False, cov_samp=50)
    print("covariance beta result", (final_cov))

    final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval, 
                                            nu_spectrum=True, cov_samp=50)
    print("covariance neutrino result", final_spect1, "\n", (final_cov1))

    """ Following steps are to calculate uncertainties """
    
    # # Calculate covariance matrix of summed best fit beta spectrum
    # covmat=convertmodel.vblist["U235"].Covariance(beta235,
    #     xval, nu_spectrum = False, samples=30)
    
    # # Calculate covariance matrix of summed converted neutrino spectrum
    # covmat_nu=convertmodel.vblist["U235"].Covariance(beta235,
    #     xval, nu_spectrum = True, samples=30)

    
    # Calculate the relative uncertainty
    relativeErr = np.zeros(len(xval))
    for i in range(len(final_spect)):
        if final_spect[i] > 0:
            relativeErr[i] = np.sqrt(final_cov[i][i])/final_spect[i]

    relativeErrNu = np.zeros(len(xval))
    for i in range(len(final_spect1)):
        if final_spect1[i] > 0:
            relativeErrNu[i] = np.sqrt(final_cov1[i][i])/final_spect1[i]

    fig = plt.figure()
    plt.ylim([-.5, .5])
    plt.fill_between(xval[0:], relativeErr[0:], -relativeErr[0:], label='beta',
                    alpha = 0.4)
    plt.fill_between(xval[0:], relativeErrNu[0:], -relativeErrNu[0:], label='neutrino',
                    alpha = 0.4)
    plt.legend()
    fig.savefig("U235_errs_new.png")

    print("beta integral", sum(final_spect))
    print("neu integral", sum(final_spect1))
    
    # generate the neutrino spectrum
    betaSpectraDB = BetaEngine(xbins=xval)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])
    
    # The following code is to test the consistency with the synthetic beta and
    # neutrino spectrum
    U235 = FissionIstp(92, 235, Ei=0)
    U235.LoadFissionDB(DB='JEFF')
    U235.LoadCorrelation(DB='JEFF')
    U235.CalcBetaSpectra(betaSpectraDB)

    # Generate the summation result
    sum1 = SumEngine(betaSpectraDB)
    sum1.AddFissionIstp(U235, "U235", 1, 0)
    # sum all spect
    sum1.CalcReactorSpectrum()

    summed_spect = sum1.spectrum


    """ Following steps are to plot the best fit spectra and compare them to 
        original spectra
    """
    
    # newxval = np.arange(2.125, 8.375, 0.25)
    # newyval = Rebin(xval, final_spect, newxval)
    # newyval1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(newxval, nu_spectrum=True, cov_samp=5)

    fig = plt.figure()
    # plt.yscale('log')
    plt.errorbar(convertmodel.betadata["U235"].x,
        convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
        label='synthetic beta')
    plt.plot(xval, final_spect, label="best fit beta")
    plt.plot(xval, summed_spect, label='synthetic neutrino')
    plt.plot(xval, final_spect1, label='best fit neutrino')
    plt.xlabel("Energy (MeV)")
    plt.ylabel("decay/fission/MeV")

    plt.legend()
    fig.savefig("bestfit_spectra.png")

    # calculating summed spectrum synthetic data as comparison
    betaspect = np.interp(xval, convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y)
    diff = (final_spect-betaspect)/betaspect
    fig = plt.figure()
    plt.ylim([-0.1, 0.1])
    plt.xlim([2, 9])
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Relative difference")
    plt.plot(xval, diff)
    plt.show()
    fig.savefig("bestfit_beta_compare.png")
    
    # Finally, the rest of this is plotting the relative errors for both the betas and the neutrinos for our given
    # Energy range
    relativeErr = final_unc/final_spect
    print(final_unc)
    relativeErrNu = final_unc1/final_spect1
    print(relativeErr, relativeErrNu)
    fig = plt.figure()
    plt.ylim([-.3, .3])
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Relative uncertainty")
    plt.fill_between(xval, -relativeErr, relativeErr, label='beta', alpha=0.4)
    plt.fill_between(xval, -relativeErrNu, relativeErrNu, label='neutrino', alpha=0.4)
    plt.legend()
    name = "U235_errors.png"
    fig.savefig(name)

    """ Following steps are to compare converted neutrino spectrum to the 
        synthetic neutrino data
    """
    diff1 = (final_spect1-summed_spect)/summed_spect*100
    fig = plt.figure()
    plt.ylim([-10, 10])
    plt.ylabel("Relative difference (\%)")
    plt.xlim([2, 9])
    plt.xlabel("Energy (MeV)")
    plt.plot(xval, diff*100, label="best-fit beta")
    plt.plot(xval, diff1, label="converted neutrino")
    plt.legend()
    fig.savefig("huber_compare.png")

    # file = open("U235_synth_compare.csv", "w")
    # file.write("E,Ne,dNe")
    # file.write("\n")
    # for i in range(len(summed_spect)):
    #     file.write(str(xval[i] * 1000) + " , " + str(summed_spect[i]) + ",0")
    # #
    # #
    # file.close()

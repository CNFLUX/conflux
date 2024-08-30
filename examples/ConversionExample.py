# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np
import os

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
    beta235 = BetaData(os.environ["CONFLUX_DB"]+"/conversionDB/U_235_e_2014.csv")
    # beta2351 = BetaData("./data/conversionDB/Synthetic_235_beta.csv")
    beta235s = BetaData(CONFLUX_DB+"/example_models//U235_synth_data_1.5_9.6.csv")
    beta239 = BetaData(os.environ["CONFLUX_DB"]+"/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData(os.environ["CONFLUX_DB"]+"/conversionDB/Pu_241_e_2014.csv")

    # Define isotopic fission yield DB to calculate average atom numbers of
    # virtual branches
    U235 = FissionIstp(92, 235, Ei=0)
    Pu239 = FissionIstp(94, 239, Ei =0)
    Pu241 = FissionIstp(94, 241, Ei = 0)

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
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)
    # convertmodel.AddBetaData(beta239, Pu239, "Pu239", 1.0)
    # convertmodel.AddBetaData(beta241, Pu241, "Pu241", 1.0)
    # Do virtual branch fitting with the defined virtual branch energy range
    xval = np.arange(0,10, 0.01)
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

    # Draw the spectra of all vertual branches
    for i in range(0, 40):
        if (sum(convertmodel.vblist["U235"].SumBranches(xval,
                thresh=i*branch_slice, nu_spectrum = False))<=0):
            continue
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
            thresh=1*branch_slice, nu_spectrum = False), fmt='--')
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval,
            thresh=1*branch_slice, nu_spectrum = True), fmt='--')
    fig.savefig("U235_convert_test.png")
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
    #

    final_spect, final_unc, final_cov = convertmodel.SummedSpectrum(xval, nu_spectrum=False, cov_samp=5)
    final_spect1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(xval, nu_spectrum=True, cov_samp=5)

    # print(final_spect)
    # print(final_spect)
    # Calculate covariance matrix of summed best fit beta spectrum
    covmat=convertmodel.vblist["U235"].Covariance(beta235,
        xval, nu_spectrum = False, samples=20)
    # Calculate covariance matrix of summed converted neutrino spectrum
    covmat_nu=convertmodel.vblist["U235"].Covariance(beta235,
        xval, nu_spectrum = True, samples=20)

    print("covariance beta result", final_cov)
    print("covariance neutrino result", final_cov1)
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
    plt.fill_between(xval, relativeErr, -relativeErr, label='beta',
                    alpha = 0.4)
    plt.fill_between(xval, relativeErrNu, -relativeErrNu, label='neutrino',
                    alpha = 0.4)
    plt.legend()
    fig.savefig("U235_errs_new.png")
    # print(final_unc)
    # print(final_cov)
    print("beta integral", sum(final_spect))
    print("neu integral",sum(final_spect1))


    newxval = np.arange(2.125, 8.375, 0.25)
    newyval = Rebin(xval, final_spect, newxval)
    newyval1, final_unc1, final_cov1 = convertmodel.SummedSpectrum(newxval, nu_spectrum=True, cov_samp=5)
    # for i in convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True):
    #     print(i)

    # testxval = np.linspace(0,200,201)
    # testyval = 1*testxval+2
    # testxout = np.linspace(0,200,21)
    # testyout = Rebin(testxval, testyval, testxout)
    # for a in (newyval1):
    #     print(a)

    fig = plt.figure()
    # plt.yscale('log')
    plt.errorbar(convertmodel.betadata["U235"].x,
        convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr,
        label='beta data')
    plt.plot(newxval, newyval1, label='neutrino rebined')
    plt.plot(xval, final_spect1, label='best fit neutrino')
    plt.legend()
    plt.show()
    fig.savefig("bestfit_spectra.png")

    xbins = np.arange(0, 8.25, 0.1)


    betaspect = np.interp(xval, convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y)
    diff = (final_spect-betaspect)/betaspect
    fig = plt.figure()
    plt.ylim([-0.1, 0.1])
    plt.xlim([2, 9])
    plt.plot(xval, diff)
    plt.show()
    fig.savefig("bestfit_beta_compare.png")

    # The following code is to test the consistency with the synthetic beta and
    # neutrino spectrum
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB(DB='JEFF')

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0, fraction=1)

    # Generate the summation result
    sum1 = SumEngine(xbins = xval)
    sum1.AddModel(model)

    # generate the neutrino spectrum
    betaSpectraDB = BetaEngine(sum1.FPYlist.keys(), xbins=xval)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    betaDBBase = BetaEngine()
    count = 0
    newlist = []

    sum1.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False)
    summed_spect = sum1.spectrum

    file = open("U235_synth_compare.csv", "w")
    file.write("E,Ne,dNe")
    file.write("\n")
    for i in range(len(summed_spect)):
        file.write(str(xval[i] * 1000) + " , " + str(summed_spect[i]) + ",0", file=file)
    #
    #
    file.close()

    #
    diff = (final_spect1-summed_spect)/summed_spect
    fig = plt.figure()
    plt.ylim([-0.3, 0.3])
    plt.xlim([2, 9])
    plt.plot(xval, diff)
    plt.show()
    fig.savefig("synthetic_compare.png")

from conflux.FPYEngine import FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData, VirtualBranch
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    betaU235 = BetaData("../data/conversionDB/U_235_e_2014.csv")
    U235 = FissionIstp(92,235)
    U235.LoadFissionDB()
    xval = np.arange(0.,10.0, 0.05)

    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(betaU235, U235, "U235", 1.0)
    convertmodel.VBfitbeta("U235", slicesize=1.0)
    fig = plt.figure()
#    plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.ylabel("(electron/neutrino)/fission/MeV")

    for i in range(0, 20):
        if not sum(convertmodel.vblist["U235"].SumBranches(xval, thresh =0.5, nu_spectrum = False))>0:
            continue
        plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = False), fmt='--')
    #Plot out the raw beta data


    # plt.errorbar(convertmodel.betadata["U235"].x, convertmodel.betadata["U235"].y, convertmodel.betadata["U235"].yerr, label='beta data')

    nu_spectrum = sum(convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = True))
    #Plot out the calculated beta spectrum
    #plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = False), label='beta')
    #Plot out the calculated neutrino spectrum
    plt.errorbar(xval, convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum = True), label='neutrino')
    plt.legend()
    plt.yscale("log")
    fig.savefig("U235_conversion_nonlog_1.png")

    for i in range(0, 20):
        if not sum(convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = True))>0:
            continue
        nu_spectrum -= convertmodel.vblist["U235"].SumBranches(xval, thresh =i*0.5, nu_spectrum = True)

    print(nu_spectrum)
    plt.clf()
    plt.plot(xval, -nu_spectrum)
    plt.yscale("log")
    plt.savefig("dumb_plot_conversion_subtracted.png")
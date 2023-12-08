from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
import numpy as np
from conflux.ConversionEngine import ConversionEngine, BetaData

if __name__ == "__main__":

    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB()
#    U235.LoadCorrelation()

    print("-----Fission Isotope Created, Covariance matrix Loaded-----")

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei=0, fraction=1)

    sumE = SumEngine(xbins=np.arange(0.0,8.0,0.1))
    sumE.AddModel(model)

    betaFDB = BetaEngine(sumE.FPYlist.keys(), xbins=np.arange(0.0,8.0,0.1))
    betaFDB.CalcBetaSpectra(nu_spectrum=False, branchErange=[0.0, 8.0])

    print("-----Summation Engine Created, Beta Spectrum Calculated-----")


    totalSpec = np.zeros(len(np.arange(0.0,8.0,0.1)))
  #--Generate the total spectrum, save it to a file

    for i in betaFDB.istplist:
        temp = betaFDB.istplist[i]
        for j in range(len(temp.spectrum)):
            if (totalSpec[j] < temp.spectrum[j]):
               totalSpec[j] = temp.spectrum[j]

    xbin = np.arange(0.0,8.0,0.1)

    file = open("U235_synth_data_2_8.csv", "w")
    file.write("E,Ne,dNe")
    file.write("\n")
    for i in range(len(totalSpec)):
        print(str(xbin[i] * 1000) + " , " + str(totalSpec[i]) + ",0", file=file)
    #
    #
    file.close()
    print("-----Synthetic data generated and saved in file-----")


    betaFDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 8.0])
    sumE.CalcReactorSpectrum(betaFDB)

    print("-----Neutrino spectrum calculated-----")


#     print(sumE.spectrum)

    #Run the conversion engine calculation
    percent = np.zeros(len(sumE.xbins))
    #Load a certain set of beta data
#    beta235 = BetaData("./U_235_e_2014.csv")
    beta235 = BetaData("./U235_synth_data_2_8.csv")
    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(beta235, U235, "U235", 1.0)
    ConvertModel.VBfitbeta(istp="U235", slicesize=0.25)
    x = sumE.xbins
#    x = np.arange(1.5,9.6,0.05)
    convertY, convertUnc, convertCov = ConvertModel.SummedSpectrum(x, cov_samp=10)


    print("-----Conversion engine run, data generated-----")

    sumY = sumE.spectrum
    sumUnc = sumE.uncertainty
    print(convertY)
    print(sumY)


    for i in range(len(percent)):
        diff = convertY[i] - sumY[i]
        percent[i] = diff/sumY[i]

    print("-----percent diff calculated-----")

    print(percent)

    #This is for plotting. Stuff above needs to work first
#     plt.plot(x, sumY, label="Summation")
#     plt.plot(x, convertY, label= "conversion")
#     plt.xlabel("Energy (in MeV)")
#     plt.ylabel("Neutrino Spectrum (neutrinos/MeV/Fission)")
#     plt.title("Slice_Size=0.25")
#     plt.legend()
#     plt.savefig("Comparison3_moreSlices4.png")
#     plt.close()
# #    plt.rcParams['text.usetex'] = True
#     plt.plot(x, percent)
#     plt.xlabel("Energy (in MeV)")
#     plt.ylabel("thing")
#     plt.title("Slice_Size=0.25")
# #    plt.plot(x, 0.2)
#     plt.savefig("Con_Sum_diff_moreSlices4.png")
#     plt.close()
#     # plt.plot(x[:-40],percent[:-40])
#     # plt.title("Slice_Size=0.15")
#     # plt.xlabel("Energy (in MeV)")
#     # plt.ylabel("fractional Diff")
#     # plt.savefig("Con_Sum_diff_moreSlices4_zoomed.png")
#

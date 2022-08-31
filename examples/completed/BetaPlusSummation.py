from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Load data into the simulation (Change this directory location to the location
    #Of U_235_e_2014.csv on the host machine)
    beta235 = BetaData("../data/conversionDB/U_235_e_2014.csv")
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB()
    U235.LoadCorrelation()
  
    #This is the conversion calculation

    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(beta235, U235, "U235", 1.0)
    ConvertModel.VBfit()

    #This is the summation calculation

    FisModel = FissionModel()
    FisModel.AddContribution(isotope=U235, Ei=0, fraction=1.0)

    result = SumEngine(spectRange=[0.0,20.0])
    result.AddModel(FisModel)

    #This is the lower range of energy for SumEngine
    betaSpectralow = BetaEngine(result.FPYlist.keys())
    #Note that I've defined the branch range to have a Q value between 0 and 1.8
    betaSpectralow.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 1.8])
    result.CalcReactorSpectrum(betaSpectralow)
    lowY = result.spectrum

    for i in result.istplist:
        print(i)

    #This is to clear out the SumEngine
    result.Clear()
    result.spectrum = np.zeros(result.bins)
    result.uncertainty = np.zeros(result.bins)

    #This is the high range of energy for SumEngine
    result.AddModel(FisModel)
    betaSpectraHigh = BetaEngine(result.FPYlist.keys())
    #Note that I've defined the branch range to have a Q value between 8 and 20
    # betaSpectraHigh.CalcBetaSpectra(nu_spectrum=True, branchErange=[8.0, 20.0])
    # result.CalcReactorSpectrum(betaSpectraHigh)
    # highY = result.spectrum

    x = np.linspace(0., 20., 200)
    convertY = ConvertModel.vblist["U235"].SumBranches(x, nu_spectrum=True)

    #Combine the spectra together
    y = []
    for i in range(len(convertY)):
        y.append(convertY[i] + lowY[i])


    fig = plt.figure()
    SumX = result.bins
    SumY = result.spectrum

    plt.plot(x,y, label = "total combined spectrum")
    #plt.plot(x, highY, label = "SumEngine Q>8.0")
    plt.plot(x, lowY, label = "SumEgine Q<1.8")
    plt.plot(x, convertY, label = "conversion Mode")
    plt.yscale('log')
    plt.title("U235 Conversion Calculation with Summation added")
    plt.legend()
    fig.savefig("HybridCalc_U235.png")

    print("Everything Worked!")

    # #Output relevant data to a text file.
    # f = open("Hybrid_Model.txt", "a")
    # f.write ("Energy, Data" + "\n")
    # for i in range(len(result.spectrum)):
    #     k = str(x[i]) + "," + str(y[i])   + "\n"
    #     f.write(k)
    # f.close()
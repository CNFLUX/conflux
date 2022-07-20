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

    #This is the summation calculation

    FisModel = FissionModel()
    FisModel.AddContribution(isotope=U235, Ei=0, fraction=1.0)

    SumEngine = SumEngine()
    SumEngine.AddModel(FisModel)

    BetaSpectraDB = BetaEngine(SumEngine.FPYlist.keys())
    BetaSpectraDB.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1)

    SumEngine.CalcReactorSpectrum(BetaSpectraDB)

    #This is the conversion calculation

    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(beta235, U235, "U235", 1.0)
    ConvertModel.VBfit()


    #This is plotting and applying cuts to the data.

    fig = plt.figure()
    convertX = np.linspace(0., 10., 200)
    convertY = ConvertModel.vblist["U235"].SumBranches(convertX, nu_spectrum = True)
    SumX = SumEngine.bins
    SumY = SumEngine.reactorSpectrum
    #cut everything below 1.8 and above 8~ MeV
    convertX= convertX[36:]
    convertY= convertY[36:]
    convertX=convertX[:124]
    convertY=convertY[:124]

    TotalX = []
    TotalY = []
    #get the summation calculation below 1.8 and above 8MeV
    TopX = SumX[80:]
    TopY = SumY[80:]
    BotX = SumX[:18]
    BotY = SumY[:18]


    #Combining all the data we've pulled from summation and conversion
    for i in range(len(BotY)):
        TotalX.append(BotX[i])
        TotalY.append(BotY[i])
    for i in range(len(convertX)):
        TotalX.append(convertX[i])
        TotalY.append(convertY[i])
    for i in range(len(TopX)):
        TotalX.append(TopX[i])
        TotalY.append(TopY[i])


    #Plot the data
    plt.plot(convertX, convertY, label="conversion mode")
    plt.plot(SumX, SumY, label="summation mode")
    plt.plot(TotalX, TotalY, "g--", label="Hybrid")
    plt.yscale('log')
    plt.legend()
    fig.savefig("Hybrid_Model.png")
    plt.close()

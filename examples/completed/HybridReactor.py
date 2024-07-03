from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Load data into the simulation (Change this directory location to the location
    #Of U_235_e_2014.csv on the host machine)
    #This is the conversion Data
    beta235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")
    beta239 = BetaData(CONFLUX_DB+"/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData(CONFLUX_DB+"/conversionDB/Pu_241_e_2014.csv")
    #This is the loaded Fission Data
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB()
    Pu239 = FissionIstp(94, 239)
    Pu239.LoadFissionDB()
    Pu241 = FissionIstp(94,241)
    Pu241.LoadFissionDB()
    U238 = FissionIstp(92,238)
    U238.LoadFissionDB()

    #This is the conversion calculation
    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(beta235, U235, "U235", .6304)
    ConvertModel.AddBetaData(beta239, Pu239, "Pu239", .2525)
    ConvertModel.AddBetaData(beta241, Pu241, "Pu241", 0.0417)
    ConvertModel.VBfitbeta(istp="U235")
    ConvertModel.VBfitbeta(istp="Pu239")
    ConvertModel.VBfitbeta(istp="Pu241")

    print("Beta Conversion calculations done")

    #This is the summation calculation
    FisModel = FissionModel()
    FisModel.AddContribution(isotope=U238, Ei=0.5, fraction=0.0754)
    result = SumEngine()
    result.AddModel(FisModel)



    #This is the U238 Summation Calculation
    betaSpectraDB = BetaEngine(result.FPYlist.keys())
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])
    result.CalcReactorSpectrum(betaSpectraDB)

    print("Summation Calculation done")

    SumY = result.spectrum

    x = np.linspace(0., 20., 200)
    convertY, uncertainty, covariance = ConvertModel.SummedSpectrum(x)

    #Combine the spectra together
    y = []
    for i in range(len(convertY)):

        y.append(convertY[i] + SumY[i])


    fig = plt.figure()
    SumY = result.spectrum

    plt.plot(x,y, label = "total combined spectrum")
    plt.plot(x, convertY, label = "Conversion Spectrum")
    plt.plot(x, SumY, label = "U238 Sum Spectrum")
    plt.yscale('log')
    plt.title("U238 Added to Conversion model for U235/Pu239/Pu241")
    plt.legend()
    fig.savefig("U238_Added.png")

    print("Everything Worked!")

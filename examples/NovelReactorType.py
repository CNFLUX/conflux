from conflux.FPYEngine import FissionIstp, FissionModel
from conflux.ConversionEngine import ConversionEngine, BetaData, VirtualBranch
from conflux.SumEngine import SumEngine
from conflux.BetaEngine import BetaEngine
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    #Load all the data into the simulation
    beta235 = BetaData("../data/conversionDB/U_235_e_2014.csv")
    beta241 = BetaData("../data/conversionDB/Pu_241_e_2014.csv")
    beta239 = BetaData("../data/conversionDB/Pu_239_e_2014.csv")
    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)
    U235.LoadFissionDB()
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()


    #Summation culation
    FisModel = FissionModel()
    FisModel.AddContribution(isotope=U235, Ei=0, fraction=0.25)
    FisModel.AddContribution(isotope=Pu239, Ei=0, fraction=0.65)
    FisModel.AddContribution(isotope=Pu241, Ei=0, fraction=0.1)

    SumEngine = SumEngine()
    SumEngine.AddModel(FisModel)
    BetaSpectraDB = BetaEngine(SumEngine.FPYlist.keys())
    BetaSpectraDB.CalcBetaSpectra(nu_spectrum=True)

    SumEngine.CalcReactorSpectrum(BetaSpectraDB)

    #Conversion Calculation
    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(beta235, U235, "U235", frac = 0.25)
    ConvertModel.AddBetaData(beta239, Pu239, "Pu239", frac = 0.65)
    ConvertModel.AddBetaData(beta241, Pu241, "Pu241", frac = 0.10)
    ConvertModel.VBfitbeta("U235")
    ConvertModel.VBfitbeta("Pu239")
    ConvertModel.VBfitbeta("Pu241")

    #Plotting

    fig = plt.figure()
    convertX = np.linspace(0., 10., 200)
    totalY, uncertainty, covariance = ConvertModel.SummedSpectrum(convertX)
    # convert235Y = ConvertModel.vblist["U235"].SumBranches(convertX, nu_spectrum = True)
    # convert239Y = ConvertModel.vblist["Pu239"].SumBranches(convertX, nu_spectrum= True)
    # convert241Y = ConvertModel.vblist["Pu241"].SumBranches(convertX, nu_spectrum=True)
    SumX = SumEngine.xbins
    SumY = SumEngine.spectrum

    # plt.plot(convertX, convert235Y, label="conversion mode 235")
    # plt.plot(convertX, convert239Y, label="conversion mode 239")
    # plt.plot(convertX, convert241Y, label="conversion mode 241")

    plt.plot(SumX, totalY, label="Conversion")
    plt.plot(SumX, SumY, label="Summation")
    plt.yscale('log')
    plt.legend()
    plt.title("Novel Reactor Type (0.25 U235, 0.65 Pu239, 0.1 Pu241)")
    fig.savefig("Novel_Reactor.png")
    plt.close()


    print("This worked")

    # 65 % Pu239,0.25 U235, 0.10 Pu241

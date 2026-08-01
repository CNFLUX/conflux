from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__=="__main__":

    #Load The Isotopic info
    x = np.arange(0., 20., 0.1)
    U235D = BetaData("../../data/conversionDB/U_235_e_2014.csv")
    Pu239D = BetaData("../../data/conversionDB/Pu_239_e_2014.csv")
    Pu241D = BetaData("../../data/conversionDB/Pu_241_e_2014.csv")
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB()
    U235.LoadCorrelation()
    Pu239 = FissionIstp(94, 239)
    Pu239.LoadFissionDB()
    Pu239.LoadCorrelation()
    Pu241 = FissionIstp(94, 241)
    Pu241.LoadFissionDB()
    Pu241.LoadCorrelation()

    U238 = FissionIstp(92, 238)
    U238.LoadFissionDB()
    U238.LoadCorrelation()


    #Load the Fission Model, add the U238 contribution to it.
    FisModel = FissionModel()
    FisModel.AddContribution(U238, 0.5, 0.25,  0.0)

    #Load a summation engine and add the fission model to it.
    Sum = SumEngine(xbins = x)
    Sum.AddModel(FisModel)

    #Load a Beta Engine and calculate the beta spectrum of the branches
    Beta = BetaEngine(Sum.FPYlist.keys(), xbins=x)
    Beta.CalcBetaSpectra(nu_spectrum=True, branchErange=[0., 20.0])

    #Calculate the total reactor spectrum
    Sum.CalcReactorSpectrum(Beta)

    #Create a conversion engine, load up the other three isotopes into it and
    #Fit the virtual branches to each isotopic spectrum
    convert = ConversionEngine()
    convert.AddBetaData(Pu241D, Pu241, "Pu241", 0.25)
    convert.AddBetaData(Pu239D, Pu239, "Pu239", 0.25)
    convert.AddBetaData(U235D, U235, "U235", 0.25)
    convert.VBfitbeta("Pu241", slicesize=0.25)
    convert.VBfitbeta("Pu239", slicesize=0.25)
    convert.VBfitbeta("U235", slicesize=0.25)

    #Pull the spectrum, uncertainty, and covariance matrix from the conversion engine.
    spectrum, uncertainty, cov = convert.SummedSpectrum(x, cov_samp = 20)


    #Plotting
    plt.plot(x, spectrum, "--", label="Converted U235/Pu239/Pu241 Spectra")
    plt.plot(x, Sum.spectrum, "--", label="Summed U238 Spectra")
    plt.plot(x, spectrum + Sum.spectrum, label="Summed Spectra")
    plt.legend()
    plt.yscale("log")
    plt.savefig("This.png")

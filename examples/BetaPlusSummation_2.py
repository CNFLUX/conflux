from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    e = np.arange(0., 20., 0.1)


    beta235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")
    U235 = FissionIstp(92, 235, 0)
    betaSpectraDB = BetaEngine(xbins = e)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])
    U235.CalcBetaSpectra(betaSpectraDB)
    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(beta235, U235, "U235", 1)
    ConvertModel.VBfitbeta(istp="U235", slice = 0.25)
    ConvSpec, ConvUnc, ConvCov = ConvertModel.SummedSpectrum(e)
    SummationEngine = SumEngine(betaSpectraDB)
    SummationEngine.AddFissionIstp(U235, "U235", count = 1)
    SummationEngine.NormalizeFP()
    SummationEngine.CalcReactorSpectrum()


    totalSpec = np.zeros(len(e))
    totalUnc = np.zeros(len(e))
    for i in range(0, 19):
        totalSpec[i] = SummationEngine.spectrum[i]
        totalUnc[i] = SummationEngine.uncertainty[i]
    for i in range(19, 76):
        totalSpec[i] = ConvSpec[i] 
        totalUnc[i]  = ConvUnc[i]
    for i in range(76, len(e)):
        totalSpec[i] = SummationEngine.spectrum[i]
        totalUnc[i] = SummationEngine.uncertainty[i]
        


    fig = plt.plot()    
    plt.errorbar(e, SummationEngine.spectrum, SummationEngine.uncertainty, fmt="-.", color = "blue", label = "Summation")
    plt.errorbar(e, ConvSpec, ConvUnc, fmt=".g", label="Conversion")
    plt.errorbar(e, totalSpec, totalUnc, fmt="-r",label="Combined")
    plt.yscale("log")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("neutrino/MeV/Fission)$")
    plt.ylim([1e-6, 10])
    plt.legend()

    # Define the x-coordinates for the vertical lines
    line_xs = [19, 76]
    
    # Add vertical lines that stop at the data points
    for line_x in line_xs:
        # Find the corresponding y-value
        y_value = totalSpec[line_x]
        plt.plot([e[line_x], e[line_x]], [0, y_value], color='k', linestyle='--')


    plt.savefig("Low_High_added.pdf")
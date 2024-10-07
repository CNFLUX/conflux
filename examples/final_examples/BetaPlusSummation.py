from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Load data into the simulation (Change this directory location to the location
    #Of U_235_e_2014.csv on the host machine, which is in the data folder)
    betaU235 = BetaData("../data/conversionDB/U_235_e_2014.csv")
    U235 = FissionIstp(92, 235, Ei = 0)
    U235.LoadFissionDB()
    U235.LoadCorrelation()
    #Initialize the BetaSpectrum, calculate the spectra of Fissile Isotopes, and 
    #Initialize the energy scale for our calculation
    e = np.arange(0, 20., 0.1)
    BetaSpectraDB = BetaEngine(xbins = e)
    BetaSpectraDB.CalcBetaSpectra(nu_spectrum=True, )

    #Next, Calculate the U235 Beta Spectrum, add it to the Summation Engine, and then
    #Calculate the total Neutrino spectrum using the summation method
    U235.CalcBetaSpectra(BetaSpectraDB)
    SummationEngine = SumEngine(BetaSpectraDB)
    SummationEngine.AddFissionIstp(U235, "U235", count = 1)
    SummationEngine.CalcReactorSpectrum()

    #We have finished the Summation Calculation, time to calculate this using the Conversion method

    #This is the conversion calculation. Initialize the Conversion model, add the U235 beta data and
    #fissile data to the model, and then fit virtual beta branches to the data.
    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(betaU235, U235, "U235", frac = 1)
    ConvertModel.VBfitbeta("U235")

    
    #Next, Calculate the Conversion Spectrum and it's uncertainties, and pull the summation spectrum
    #And uncertainties from the summation engine.
    convertSpec, convertUnc, convertCov = ConvertModel.SummedSpectrum(e, nu_spectrum=True, cov_samp=5)
    SumSpec = SummationEngine.spectrum
    SumUnc = SummationEngine.uncertainty

    #Now, we plot and apply cuts to the data. 
    fig = plt.figure()

    #I create an empty array to store my combined spectrum, and then I do the following:
    #Since the conversion calculation cannot predict the neutrino spectrum below 1.8MeV (IBD threshold)
    #We take the summation spectrum below 1.8MeV and "tack it" onto the conversion spectrum above 1.8MeV
    #And below 8MeV. I then take the summation spectrum and "tack" the spectrum above 8MeV onto my overall
    #spectrum. I do this for the uncertainties as well
    TotalY = []
    TotalUnc = []
    for i in range(0,18):
        TotalY.append(SumSpec[i])
        TotalUnc.append(SumUnc[i])
    for i in range(18, 80):
        TotalY.append(convertSpec[i])
        TotalUnc.append(convertUnc[i])
    for i in range(80, 200):
        TotalY.append(SumSpec[i])
        TotalUnc.append(SumUnc[i])



    #Finally, I Plot the data and save it as a hybrid model
    plt.errorbar(e, convertSpec, yerr = convertUnc, label="conversion mode")
    plt.errorbar(e, SumSpec, yerr = SumUnc, label="summation mode")
    plt.errorbar(e, TotalY, yerr = TotalUnc, fmt = "--", ecolor = "pink", label="Hybrid")
    plt.yscale('log')
    plt.xlabel("Energy (in MeV)")
    plt.ylabel("Spectrum (neutrinos/MeV/fission)")
    plt.legend()
    fig.savefig("Hybrid_Model.png")
    plt.close()

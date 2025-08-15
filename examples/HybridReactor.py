from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np



#Here, we will calculate the total spectrum for a Reactor using both the Summation and Conversion Engines
#The conversion engine will be used for U235, Pu239, and Pu241. Since the Beta Database for U238 is incomplete,
#We will calculate its' spectrum using the ab-initio method before adding it to the overall spectrum.
if __name__ == "__main__":
    #First, I load up the beta data. This will be in the data folder in the CONFLUX installation directory
    beta235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")
    beta239 = BetaData(CONFLUX_DB+"/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData(CONFLUX_DB+"/conversionDB/Pu_241_e_2014.csv")
    #Next, I load in U235, U238, Pu239, and Pu241 Fission Isotopes. I do not need to call LoadFissionDB(), as
    #The code will automatically load it when the FissionIstps are initialized
    U235 = FissionIstp(92, 235, 0)
    Pu239 = FissionIstp(94, 239, 0)
    Pu241 = FissionIstp(94,241, 0)
    U238 = FissionIstp(92,238, 0.5)

    #While I'm at it, I will initailize an energy range, and calculate the neutrino spectrum of my Fission Products
    #With the BetaEngine
    e = np.arange(0., 14., 0.1)
    betaSpectraDB = BetaEngine(xbins = e)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    #I will then Calculate the beta spectrum for U235, Pu239, and Pu241 and U238
    U235.CalcBetaSpectra(betaSpectraDB)
    Pu239.CalcBetaSpectra(betaSpectraDB)
    Pu241.CalcBetaSpectra(betaSpectraDB)
    U238.CalcBetaSpectra(betaSpectraDB)

    #Next, I will run the conversion calculation. I add U235, Pu239, and Pu241 to the Conversion model and then
    #Individually fit their virtual branches. Finally, I calculate the total spectrum by summing all spectra
    #That are in the conversion model together.
    ConvertModel = ConversionEngine()
    ConvertModel.AddBetaData(beta235, U235, "U235", .6304)
    ConvertModel.AddBetaData(beta239, Pu239, "Pu239", .2525)
    ConvertModel.AddBetaData(beta241, Pu241, "Pu241", 0.0417)
    ConvertModel.VBfitbeta(istp="U235")
    ConvertModel.VBfitbeta(istp="Pu239")
    ConvertModel.VBfitbeta(istp="Pu241")
    ConvSpec, ConvUnc, ConvCov = ConvertModel.SummedSpectrum(e, cov_samp=1)


    #Now that I'm done with the conversion calculation, let me focus on the summation engine
    #I will Load up a Summation model, add our U238 FissionIstp to it, and then calculate the total
    #Neutrino spectrum of U238
    SummationEngine = SumEngine(betaSpectraDB)
    SummationEngine.AddFissionIstp(U238, "U238", count = 0.0754)
    SummationEngine.CalcReactorSpectrum()

    #Now, I need to add both the spectrum from the U238 Calculation, as well as the Conversion Calculation
    #I add them in a very straightforward manner. Note, We have already given the Fractional contribution
    #Each of these isotopes has to the overall reactor spectrum.
    y = []
    y_err = []
    for i in range(len(e)):
        y.append(ConvSpec[i] + SummationEngine.spectrum[i])
        y_err.append(ConvUnc[i] + SummationEngine.uncertainty[i])

    #Finally, I go ahead and plot my individual engine spectra, as well as the total reactor spectrum, complete
    #With calculated uncertainties.
    fig = plt.figure()
    plt.errorbar(e[0:150],y[0:150], yerr = SummationEngine.uncertainty[0:150] + ConvUnc[0:150], label = "total combined spectrum")
    plt.errorbar(e[0:150], ConvSpec[0:150], yerr = ConvUnc[0:150], label = "Conversion Spectrum")
    plt.errorbar(e[0:150], SummationEngine.spectrum[0:150],yerr = SummationEngine.uncertainty[0:150], label = "U238 Sum Spectrum")
    plt.ylabel(r"${\nu}_e/MeV/Fission$")
    plt.xlabel("E (in MeV)")
    plt.legend()
    fig.savefig("Hybrid_reactor.pdf")
    plt.yscale('log')
    fig.savefig("Hybrid_reactor_log.pdf")
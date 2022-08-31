# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    #Beta Conversion Data
    beta235 = BetaData("../data/conversionDB/U_235_e_2014.csv")
    beta239 = BetaData("../data/conversionDB/Pu_239_e_2014.csv")

    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)

    U235.LoadFissionDB()
    Pu239.LoadFissionDB()

    #Run the Beta Conversion Engine
    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 0.6)
    convertmodel.AddBetaData(beta239, Pu239, "Pu239", 0.4)

    convertmodel.VBfit(0.25)


    #Cast spectrum and energy into more manageable lists
    xval = np.linspace(0., 10., 200)
    totalSpec = convertmodel.SummedSpectrum(xval)
    U235Spec = convertmodel.vblist["U235"].SumBranches(xval, nu_spectrum=True)
    Pu239Spec = convertmodel.vblist["Pu239"].SumBranches(xval, nu_spectrum=True)


    #Plotting
    fig = plt.plot()
    plt.yscale('log')
    plt.ylabel("neutrino/fission/MeV")
    plt.xlabel("E (in MeV)")
    plt.plot(xval, totalSpec, label = "60% U235, 40% Pu239")
    plt.plot(xval, U235Spec, label="U235")
    plt.plot(xval, Pu239Spec, label= "Pu239")
    plt.title("Comparison between U235, Pu239 and a 60/40 mix of U235/Pu239")
    plt.legend()

    plt.savefig("60-40")



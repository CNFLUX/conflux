# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    #Beta Conversion Data
    beta235 = BetaData("../../data/conversionDB/U_235_e_2014.csv")
    beta239 = BetaData("../../data/conversionDB/Pu_239_e_2014.csv")

    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)

    U235.LoadFissionDB()
    Pu239.LoadFissionDB()

    #Run the Beta Conversion Engine
    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 0.6)
    convertmodel.AddBetaData(beta239, Pu239, "Pu239", 0.4)

    convertmodel.VBfitbeta("U235", 0.25)
    convertmodel.VBfitbeta("Pu239", 0.25)

    othermodel = ConversionEngine()
    othermodel.AddBetaData(beta235, U235, "U235", 1.0)
    othermodel.VBfitbeta("U235", 0.25)



    #Cast spectrum and energy into more manageable lists
    xval = np.linspace(0., 10., 200)


    #Something is wrong here, fix it~
    totalSpec, uncertainty, covariance = convertmodel.SummedSpectrum(xval)
    otherSpec = othermodel.vblist["U235"].SumBranches(xval, nu_spectrum=True)

    # for i in otherSpec:
    #     print(i)

    #Fractional difference
    comp = []
    for i in range(len(totalSpec)):
        diff = totalSpec[i] - otherSpec[i]
        avg = (totalSpec[i] + otherSpec[i])/2.
        print(avg)
        tot = diff/avg
        # if avg == 0:
        #     tot = 0
        comp.append(tot)


    # #Plotting
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.set_yscale('log')
    ax1.set(ylabel = "(electron/neutrino)/fission/MeV")
    ax1.plot(xval, totalSpec, label = "60% U235, 40% Pu239")
    ax1.plot(xval, otherSpec, label="100% U235")
    ax2.plot(xval, comp)
    ax2.set(xlabel = "E (MeV)", ylabel = "Fractional difference")
    ax1.set_title("Comparison between pure U235 and a 60/40 mix of U235/Pu239")
    ax1.legend()

    plt.savefig("60-40")



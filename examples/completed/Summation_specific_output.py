from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import os
import matplotlib.pyplot as plt
import numpy as np



if __name__ == "__main__":

    #Initiate the Beta Isotope for the Jeff Library
    U235J = FissionIstp(92, 235)
    U235J.LoadFissionDB(defaultDB="JEFF")
    #Initialize the fission Model, and add the U235 contribution to it
    modelJ = FissionModel()
    modelJ.AddContribution(isotope=U235J, Ei=0, fraction=1)

    #Initialize the JEFF Summation calculation
    Jeff = SumEngine()
    Jeff.AddModel(modelJ)
    #Run the Summation calculation for JEFF
    betaSpectraDBJ = BetaEngine(Jeff.FPYlist.keys())
    betaSpectraDBJ.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    Jeff.CalcReactorSpectrum(betaSpectraDBJ)

    #Initialize the Beta Isotope for the ENDF Library
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB(defaultDB="ENDF")
    #Initialize the fission Model, and add the U235 contribution to it
    model = FissionModel()
    model.AddContribution(isotope=U235, Ei=0, fraction=1)
    #Initialize the ENDF Summation Calculation
    Endf = SumEngine()
    Endf.AddModel(model)
    # Run the Summation calculation for ENDF
    betaSpectraDB = BetaEngine(Endf.FPYlist.keys())
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    Endf.CalcReactorSpectrum(betaSpectraDB)

    #Calculate the fractional difference between the
    #used summation databases.
    summedDiff = []
    for i in range(len(Jeff.spectrum)):
        added =  Endf.spectrum[i] - Jeff.spectrum[i]
        average = (Jeff.spectrum[i] + Endf.spectrum[i])/2.
        total = added/average
        summedDiff.append(total)


    #Plotting

    x = np.linspace(0.,20.,200)

    fig, (ax1, ax2) = plt.subplots(2)
    ax1.plot(x, Jeff.spectrum, label="JEFF")
    ax1.plot(x, Endf.spectrum, label='ENDF')
    ax1.set_yscale("log")
    ax1.set_xlim(0,14.0)
    ax2.set(xlabel = "E (in MeV)", ylabel = "Fractional Difference")
    ax1.set(ylabel = "(electrons/neutrinos) /fission/MeV")

    ax1.legend()
    ax2.plot(x, summedDiff)
    ax1.set_title("JEFF vs ENDF for U235")
    plt.savefig("JEFFvENDF")

    print("This worked")
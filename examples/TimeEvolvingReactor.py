from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    # Load in the Fission Data for our fissile isotopes.
    # Also load up the energy range array
    
    df = pd.read_csv(CONFLUX_DB+'/example_models/timeEvolvingReactor.csv')
    
    U235 = FissionIstp(92, 235, Ei=0)
    U235.LoadFissionDB()
    U238 = FissionIstp(92, 238, Ei=0.5)
    U238.LoadFissionDB()
    Pu239 = FissionIstp(94, 239, Ei=0)
    Pu239.LoadFissionDB()
    Pu241 = FissionIstp(94, 241, Ei=0)
    Pu241.LoadFissionDB()

    e = np.arange(0,15., 0.1)
    binwidth = e[1]-e[0]

    #Load up the BetaSpectrum, and Calculate the beta spectrum for our fission products.
    betaSpectraDB = BetaEngine(xbins=e)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0,15])
    U235.CalcBetaSpectra(betaSpectraDB)
    U238.CalcBetaSpectra(betaSpectraDB)
    Pu239.CalcBetaSpectra(betaSpectraDB)
    Pu241.CalcBetaSpectra(betaSpectraDB)

    #Create a summation engine, load in the fissile isotopes with an initial contribution of 0
    SummationEngine = SumEngine(betaSpectraDB)
    SummationEngine.AddFissionIstp(U235, "U235", count = 0)
    SummationEngine.AddFissionIstp(U238, "U238", count = 0)
    SummationEngine.AddFissionIstp(Pu239, "Pu239", count = 0)
    SummationEngine.AddFissionIstp(Pu241, "Pu241", count = 0)

    #Next, Open up the reactor output file that contains the % - fissile contribution
    df = pd.read_csv(CONFLUX_DB+'/example_models/timeEvolvingReactor.csv')

    #read out the first line which should be the isotope names that are being fed into our simulation
    #As well as the # of Days (This is the Header)

    fig, (ax1, ax2) = plt.subplots(2)

    days = []
    ratio_to_day_0 = []

    #Next, I iterate through the lines in the csv file
    for index, row in df.iterrows():
        #Figure out what day we are doing the calculation for
        name = (row["Days"])
        # We do simulation calculations for only specific days
        if name not in [0.1, 20, 300, 800, 1400]:
            continue

        #Edit the fractional contribution from each fissile isotope, and then calculate the 
        #Total Reactor spectrum.
        SummationEngine.EditContribution(istpname="U235", count = row["U235"], d_count = 0)
        SummationEngine.EditContribution(istpname="U238", count = row["U238"], d_count = 0)
        SummationEngine.EditContribution(istpname="Pu239", count = row["Pu239"], d_count = 0)
        SummationEngine.EditContribution(istpname="Pu241", count = row["Pu241"], d_count = 0)
        SummationEngine.CalcReactorSpectrum()

        #Plot the spectrum at the current timestep
        sumX = SummationEngine.xbins
        sumY = SummationEngine.spectrum
        days.append(name)
        ratio_to_day_0.append(sum(sumY)*binwidth)
        Labelmaker = "Reactor on for " + str(name) + " days"
        ax1.plot(sumX, sumY, label = Labelmaker)

    #Plotting
    ax1.legend()
    ax1.set(xlabel="E (MeV)", ylabel="neurtino/MeV")
    ax2.plot(days, np.array(ratio_to_day_0)/ratio_to_day_0[0]*100, 'bo')
    ax2.set_xscale('log')
    ax2.set(xlabel = "Days since reactor on", ylabel=r"${\phi}_x / {\phi}_0$ (%)" )

    plt.savefig("time_dependent_neutrino_flux.pdf")

from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

from copy import deepcopy

if __name__ == "__main__":
    # Load in the Fission Data for our fissile isotopes.
    # Also load up the energy range array

    # Use local data path if running from examples directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_data_path = os.path.join(script_dir, '../data/example_models/timeEvolvingReactor.csv')

    if os.path.exists(local_data_path):
        csv_path = local_data_path
    else:
        csv_path = CONFLUX_DB+'/example_models/timeEvolvingReactor.csv'

    df = pd.read_csv(csv_path)
    
    U235 = FissionIstp(92, 235, Ei=0)
    U235.LoadFissionDB()
    U238 = FissionIstp(92, 238, Ei=0.5)
    U238.LoadFissionDB()
    Pu239 = FissionIstp(94, 239, Ei=0)
    Pu239.LoadFissionDB()
    Pu241 = FissionIstp(94, 241, Ei=0)
    Pu241.LoadFissionDB()

    # Save FPY data before CalcBetaSpectra modifies it
    U235_FPY = deepcopy(U235.FPYlist)
    U238_FPY = deepcopy(U238.FPYlist)
    Pu239_FPY = deepcopy(Pu239.FPYlist)
    Pu241_FPY = deepcopy(Pu241.FPYlist)

    e = np.arange(0,15., 0.1)
    binwidth = e[1]-e[0]

    #Load up the BetaSpectrum, and Calculate the beta spectrum for our fission products.
    # # Calculate beta spectra of all beta unstable isotopes
    xbins = e
    betaSpectraDB = BetaEngine(xbins=xbins)
    filename = "beta_spectra.csv"
    try:
        with open(filename, "r") as file:
            betaSpectraDB.LoadFile(filename)
    except FileNotFoundError:
        print("File not found. Creating the file.")
        betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
        betaSpectraDB.SaveToFile(filename)
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
    df = pd.read_csv(csv_path)

    #read out the first line which should be the isotope names that are being fed into our simulation
    #As well as the # of Days (This is the Header)

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.set_figheight(10)
    fig.set_figwidth(8)

    days = []
    ratio_to_day_0 = []

    count = 0
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
        if count == 0:
            first_spec = deepcopy(sumY)
        days.append(name)
        ratio_to_day_0.append(sum(sumY)*binwidth)
        Labelmaker = "Reactor on for " + str(name) + " days"
        ax1.plot(sumX, sumY, label = Labelmaker)
        ax2.plot(sumX, sumY/first_spec, label=Labelmaker)
        count+=1

    #Plotting
    ax1.legend()
    ax1.set(xlabel="E (MeV)", ylabel="neurtino/MeV")

    ax2.set(xlabel="E (MeV)", ylabel="relative to first period")

    ax3.plot(days, np.array(ratio_to_day_0)/ratio_to_day_0[0]*100, 'bo')
    ax3.set_xscale('log')
    ax3.set(xlabel = "Days since reactor on", ylabel=r"${\phi}_x / {\phi}_0$ (%)" )

    plt.savefig("time_dependent_neutrino_flux.pdf")

    # Create fission product yield plots for each fissile isotope
    fig2, ((ax4, ax5), (ax6, ax7)) = plt.subplots(2, 2, figsize=(12, 10))

    # Plot fission product yields vs mass number for each isotope
    # Use saved FPY data (before CalcBetaSpectra modified it)
    fissile_isotopes = [
        (U235_FPY, "U-235", U235.Ei, ax4),
        (U238_FPY, "U-238", U238.Ei, ax5),
        (Pu239_FPY, "Pu-239", Pu239.Ei, ax6),
        (Pu241_FPY, "Pu-241", Pu241.Ei, ax7)
    ]

    for FPYlist, name, Ei, ax in fissile_isotopes:
        # Extract mass numbers and yields from FPYlist
        mass_numbers = []
        yields = []

        for FPZAI, nuclide in FPYlist.items():
            A = (FPZAI // 10) % 1000  # Extract mass number
            mass_numbers.append(A)
            yields.append(nuclide.y)  # y is the yield value

        # Sort by mass number for cleaner plotting
        sorted_data = sorted(zip(mass_numbers, yields))
        mass_numbers, yields = zip(*sorted_data) if sorted_data else ([], [])

        print(f"{name}: {len(mass_numbers)} fission products, yield range: {min(yields) if yields else 0:.2e} to {max(yields) if yields else 0:.2e}")

        # Plot as stem plot or scatter
        if yields:
            ax.semilogy(mass_numbers, yields, 'o-', markersize=3, linewidth=0.5)
            ax.set_ylim([1e-18, 1e0])  # Set consistent y-axis range
        ax.set_xlabel('Mass Number (A)')
        ax.set_ylabel('Fission Yield')
        ax.set_title(f'{name} (Ei={Ei} MeV)')
        ax.grid(True, alpha=0.3)
        ax.set_xlim([60, 180])

    plt.tight_layout()
    plt.savefig("fission_product_yields.pdf")

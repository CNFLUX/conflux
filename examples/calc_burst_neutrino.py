# This code calculate the burst fission neutrino flux from the 
# thermal fission of U235, Pu239, and Pu241, as well as the fast
# neutron fission of U238. 
# The code contains the genenration of the default neutrino spectra
# based on the ENSDF nuclear database. If the default neutrio 
# spectra can be found in $CONFLUX_DB, the calculation will be 
# skipped.
# The second part of the code contains the calculation of the burst
# fission neutrino flux from 1e-10 to 1e9 time window in the log
# scale with 100 steps. It also looks for existing calculation. 
# If the existing calculation data can be found in 
# $CONFLUX_DB/example_models, the calculation will also be skipped.
import numpy as np
import pandas as pd

# conflux modules
from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionIstp

xbins = np.arange(0, 15, 0.01)

# Calculate beta spectra of all beta unstable isotopes
betaSpectraDB = BetaEngine(xbins=xbins)
# Loading the pre-calculated neutrino spectra
filename = f"{CONFLUX_DB}/default_neutrino_spectra.csv"
try:
    with open(filename, "r") as file:
        betaSpectraDB.LoadFile(filename)
except FileNotFoundError:
    print("File not found. Creating the file.")
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
    betaSpectraDB.SaveToFile(filename)
print("File created.")

# Declare all fissile isotopes, using the JEFF database
# NOTICE: if using the ENDF, the fast neutron energy should 
# be set to 0.5 MeV
FPY_DB_name = "JEFF"
fissile_istps = {}

fissile_istps["U235_thermal"] = FissionIstp(92, 235, 0., DB=FPY_DB_name, IFPY=True)
fissile_istps["U235_thermal"].LoadFissionDB(DB = FPY_DB_name)
fissile_istps["U235_thermal"].LoadCorrelation(DB = FPY_DB_name)

fissile_istps["U238_fast"] = FissionIstp(92, 238, 0.4, DB=FPY_DB_name, IFPY=True)
fissile_istps["U238_fast"].LoadFissionDB(DB = FPY_DB_name)
fissile_istps["U238_fast"].LoadCorrelation(DB = FPY_DB_name)

fissile_istps["Pu239_thermal"] = FissionIstp(94, 239, 0., DB=FPY_DB_name, IFPY=True)
fissile_istps["Pu239_thermal"].LoadFissionDB(DB = FPY_DB_name)
fissile_istps["Pu239_thermal"].LoadCorrelation(DB = FPY_DB_name)

fissile_istps["Pu241_thermal"] = FissionIstp(94, 241, 0., DB=FPY_DB_name, IFPY=True)
fissile_istps["Pu241_thermal"].LoadFissionDB(DB = FPY_DB_name)
fissile_istps["Pu241_thermal"].LoadCorrelation(DB = FPY_DB_name)

logedges = np.logspace(-9, 9, num=100)
# logedges = np.insert(logedges, 0, 0.0)

# calculate the burst fission of every fissile isotopes
for key in fissile_istps.keys():
    time_labels = [f"{t}" for t in logedges]    
    col_labels  = [f"{j}" for j in (xbins)] 
    spectrum = np.empty((len(logedges), len(xbins)))
    uncertainty = np.empty((len(logedges), len(xbins)))
    outname = f'{key}_burst_fission_neutrino'
    spectfile = f'{CONFLUX_DB}/example_models/{outname}.csv'
    uncfile = f'{CONFLUX_DB}/example_models/{outname}_unc.csv'
    try:
        with open(spectfile, "r") as file:
            continue
    except FileNotFoundError:
        print(f"Calculting burst fission spectrum of {key}")

        for i in range(len(logedges)):
            tau = logedges[i]
            fissile_istps[key].CalcBetaSpectra(betaSpectraDB, 
                                               processMissing=False,  
                                               time=tau, silent=False)
            spectrum[i] = (fissile_istps[key].spectrum)
            uncertainty[i] = (fissile_istps[key].uncertainty)

        df_spect = pd.DataFrame(spectrum,
                        index=time_labels,
                        columns=col_labels)
        df_spect.index.name = 'time_s'   

        df_spect.to_csv(spectfile)

        df_unc = pd.DataFrame(uncertainty,
                        index=time_labels,
                        columns=col_labels)
        df_unc.index.name = 'time_s'   

        df_unc.to_csv(uncfile)

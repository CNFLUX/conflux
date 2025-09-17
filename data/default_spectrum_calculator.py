# Prepackaged program to calculate precalculated the neutrino spectra of all isotopes in ENSDF.
# Binning is set to 0 to 15 MeV with 0.01 MeV steps. 
# Beta decay parameters are set to default values.

import numpy as np

# conflux modules
from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.BetaPlusEngine import BetaPlusEngine

xbins = np.arange(0, 15, 0.01)

# Calculate beta spectra of all beta unstable isotopes
betaSpectraDB = BetaEngine(xbins=xbins)
filename = f"{CONFLUX_DB}/default_neutrino_spectra_V2.csv"
betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
betaSpectraDB.SaveToFile(filename)
print(f"File {filename} created.")

# Calculate beta spectra of all beta unstable isotopes
betaplusSpectraDB = BetaPlusEngine(xbins=xbins)
filename2 = f"{CONFLUX_DB}/default_ec_neutrino_spectra_V2.csv"
betaplusSpectraDB.CalcBetaSpectra(nu_spectrum=True)
betaplusSpectraDB.SaveToFile(filename2)
print(f"File {filename2} created.")

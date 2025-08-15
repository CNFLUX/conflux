# Prepackaged program to calculate precalculated the neutrino spectra of all isotopes in ENSDF.
# Binning is set to 0 to 15 MeV with 0.01 MeV steps. 
# Beta decay parameters are set to default values.

# conflux modules
from conflux.BetaEngine import BetaEngine, CONFLUX_DB

xbins = np.arange(0, 15, 0.01)

# Calculate beta spectra of all beta unstable isotopes
betaSpectraDB = BetaEngine(xbins=xbins)
filename = f"{CONFLUX_DB}/default_neutrino_spectra.csv"
try:
    with open(filename, "r") as file:
        betaSpectraDB.LoadFile(filename)
except FileNotFoundError:
    print("File not found. Creating the file.")
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
    betaSpectraDB.SaveToFile(filename)
print("File created.")

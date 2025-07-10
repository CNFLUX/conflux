import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# conflux modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp
from conflux.IBDXSection import ibd_xsection_cm2

E_fission = 3.2e-11                 # energy per fission, joules
E_TNT = 4.184e9                     # energy per ton TNT, joules
fissions_per_ton = E_TNT/E_fission  # number of fissions per ton TNT

# 1000 ton TNT equivalent fission burst
m_TNT = 1000
test_flux = m_TNT * fissions_per_ton

xbins = np.arange(0, 13, 0.1)

# Calculate beta spectra of all beta unstable isotopes
betaSpectraDB = BetaEngine(xbins=xbins)
filename = "fission_burst.csv"
try:
    with open(filename, "r") as file:
        betaSpectraDB.LoadFile(filename)
except FileNotFoundError:
    print("File not found. Creating the file.")
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
    betaSpectraDB.SaveToFile(filename)
print("File created.")

# Load the independent fission product yields from the JEFF database
# Ei is set to be 0.4 MeV
U235 = FissionIstp(92, 235, 0., DB='JEFF', IFPY=False)
U235.LoadFissionDB(DB = 'JEFF')
U235.LoadCorrelation(DB = 'JEFF')
U235.CalcBetaSpectraOld(betaSpectraDB)
totalspect = U235.spectrum # * test_flux # the full spectrum summed with all fission products
print(sum(U235.spectrum))

totalspect_ibd = totalspect*ibd_xsection_cm2(xbins) # the IBD detector spectrum

# define the time windows of the calculation
windows = np.linspace(0, 100, 101)
windows_log = np.logspace(-2, 4, 7)
edges = np.concatenate(([0.0], np.logspace(-2,6,9)))

# define lists of spectra in different time windows
spect_time = []
spect_sum_t = []
spect_sum_t_thresh = []
spect_sum_t_ibd = []
cumuspect = np.zeros(len(xbins))

U2351 = FissionIstp(92, 235, 0., DB='JEFF', IFPY=True)
U2351.LoadFissionDB(DB = 'JEFF')
U2351.LoadCorrelation(DB = 'JEFF')
eval_count = 0
def spectrum_t(u):
    global eval_count
    eval_count += 1
    if eval_count % 1000 == 0:
        print(f"g called {eval_count}×, latest u = {u:.6f}")
    t = 10**u
    U2351.CalcBetaSpectra(betaSpectraDB, processMissing=False,  time=t, silent=True)
    return U2351.spectrum*t*np.log(10)

# logedges = np.array([-30, -2, -1, 0, 1, 2, 3, 4, 5, 6])
logedges = np.array([-30, -1, 0, 1, 2, 3, 4, 5])
linedges = np.power(10., logedges)

fig, ax = plt.subplots()
# generate neutrino spectra in different time window after ignition
from scipy.integrate import quad_vec
for t0, t1 in tqdm(zip(logedges[:-1], logedges[1:])):
    
    # calculate the beta spectra within different windows after the fission 
    spect = quad_vec(spectrum_t, t0, t1, epsabs=1e-6, epsrel=1e-6)[0] #*test_flux
    # print("1e%g s - 1e%g s"%(t0, t1))

    # add spectrum of each window to a cumulative spectrum
    cumuspect += spect
    spect_time.append(spect)
    spect_sum_t.append(100 * sum(cumuspect) / sum(totalspect))
    spect_sum_t_thresh.append(100 * sum(cumuspect[xbins > 1.8]) / sum(totalspect[xbins > 1.8]))
    spect_sum_t_ibd.append(100 * sum(cumuspect*ibd_xsection_cm2(xbins)) / sum(totalspect_ibd))

    ax.plot(xbins, cumuspect, label="1E%g s - 1E%g s"%(t0, t1))
    
ax.set(xlabel='E (MeV)', ylabel='neutrinos/MeV from %g ton yield'%m_TNT)
ax.plot(xbins, totalspect, label="total")
# ax.set_yscale("log")
ax.legend()
fig.savefig("235U_ENDF_jeff_0.4_MeV_time_Old.pdf")

# Drawing the total neutrino flux with respect to time 
fig, ax = plt.subplots()
ax.set(xlabel='log10 time (s)', ylabel='Percent', title='U-235 neutrino flux over time')
ax.plot(linedges[1:], spect_sum_t, label='total neutrino', color='red')
ax.plot(linedges[1:], spect_sum_t_thresh, label='neutrino above 1.8 MeV', color='green')
ax.plot(linedges[1:], spect_sum_t_ibd, label='IBDs', color='blue')
xval = np.array([1, 10, 100])
yval = np.interp(xval, linedges[1:], spect_sum_t_ibd)
# Mark the data point on the curve
for j in range(len(xval)):
    ax.scatter(xval[j], yval[j], color='red')
    # Draw a line to the axes
    ax.axhline(y=yval[j], color='gray', linestyle='--', xmin=0, xmax=xval[j])
    ax.axvline(x=xval[j], color='gray', linestyle='--', ymin=0, ymax=yval[j])

plt.legend()
plt.xscale('log')
fig.savefig("235U_ENDF_jeff_0.4_MeV_ratevstime_log_Old.pdf")


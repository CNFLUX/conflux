import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd

from scipy.interpolate import interp1d, PchipInterpolator, RegularGridInterpolator
from scipy.integrate import simpson

# conflux modules
from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionIstp
from conflux.IBDXSection import ibd_xsection_cm2

binwidth = 0.01
xbins = np.arange(0, 15, binwidth)

def fission_rate_t(timex, x, y):
    R_i_interp = interp1d(
        x,
        y,
        kind='cubic',           
        fill_value='extrapolate',
        bounds_error=False
    )
    return R_i_interp(timex)

def spectrum_t(timex, seconds, spectrum):
    timeslice = interp1d(
        seconds,
        spectrum,          # shape (N_t, N_E)
        kind="cubic",
        axis=0,
        fill_value='extrapolate',
        bounds_error=False
    )
    return timeslice(timex)

def plot_integrand(timex, Eidx, tau, fission_time, fission_rate, beta_time, beta_spec):
    Rvals = fission_rate_t(tau, fission_time, fission_rate)  # (N+1,)
    Svals = spectrum_t(timex - tau, beta_time, beta_spec)    # (N+1, N_E)
    integrand = Rvals * Svals[:,Eidx]
    plt.loglog(tau, integrand)
    plt.xlabel("τ [s]")
    plt.ylabel(f"integrand at E-bin {Eidx}")
    plt.grid(True, which="both")
    plt.show()

def convolve_via_s_log(t_x, 
                       fission_time, 
                       fission_rate, 
                       beta_time, 
                       beta_spec, 
                       N=500, 
                       s_min=1e-9):
    # build a log‐spaced grid in s = t_x - tau
    u = np.linspace(np.log10(s_min), np.log10(t_x), N+1)  # uniform in log(s)
    s = 10**u                                             # let s = 10^u a log scale integration
    # trapezoid weights in u
    w_u = np.full_like(u, (u[-1] - u[0]) / N)            
    w_u[0] *= 0.5
    w_u[-1] *= 0.5
    jac_s = np.log(10) * s                                # ds = log(10)*10**u * du

    # map back to τ = t_x - s
    tau = t_x - s

    # 3) evaluate integrand
    Rvals = fission_rate_t(tau, fission_time, fission_rate)    # shape (N+1,)
    Svals = spectrum_t(s, beta_time, beta_spec)                # shape (N+1, N_E)

    integrand = Rvals[:,None] * Svals * jac_s[:,None]
    # spectrum  = np.sum(integrand * w_u[:,None], axis=0)
    spectrum = simpson(y=integrand, x=u, axis=0)
    return spectrum
    
if __name__ == "__main__":
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

    # Loading the example fission rate data
    df_fission = pd.read_csv(CONFLUX_DB+'/example_models/timeEvolvingReactor.csv')
    
    rates = {"U235":[], "U238":[], "Pu239":[], "Pu241":[]}
    for key in rates.keys():
        fission_time = []
        days = []
        for index, row in df_fission.iterrows():
            name = (row["Days"])
            days.append(name)
            fission_time.append(name*3600*24)
            rates[key].append(row[key])
        rates[key]=np.array(rates[key])
        fission_time = np.array(fission_time)

    spectra_files = {
        "U235": f"{CONFLUX_DB}/example_models/U235_thermal_burst_fission_neutrino.csv",
        "U238": f"{CONFLUX_DB}/example_models/U238_fast_burst_fission_neutrino.csv",
        "Pu239": f"{CONFLUX_DB}/example_models/Pu239_thermal_burst_fission_neutrino.csv",
        "Pu241": f"{CONFLUX_DB}/example_models/Pu241_thermal_burst_fission_neutrino.csv",
    }
    spectrum_E_t = {}
    totalflux = {}
    spectrum_day = {}
    summed_flux = []

    days = np.linspace(0, 100)
    plt.figure()
    for key, item in spectra_files.items():
        df_spectrum = pd.read_csv(item)
        beta_time = df_spectrum.iloc[:, 0].values
        energies = df_spectrum.columns[1:].astype(float).values
        spectrum_E_t[key] = df_spectrum.iloc[:, 1:].values
 
        totalflux[key] = []
        rates_1 = np.zeros(len(rates[key]))+1
        # fig = plt.figure()
        spectrum_day[key] = []
        for day in days:
            timex = day*3600*24
            spectrum = convolve_via_s_log(timex,                                     
                                        fission_time=fission_time, 
                                        fission_rate=rates[key], 
                                        beta_time=beta_time, 
                                        beta_spec=spectrum_E_t[key], 
                                        N=1000, 
                                        )
            spectrum_day[key].append(spectrum)
            totalflux[key].append(sum(spectrum[180:])*0.01)

        plt.plot(days, totalflux[key], label=key, linestyle="-" )
 

    totalspect = []
    for i in range(len(days)):
        spectrum = np.zeros(len(xbins))
        for key, item in spectra_files.items():
            spectrum += spectrum_day[key][i]
        totalspect.append(sum(spectrum)*0.01)
    plt.plot(days, totalspect, label="total")
    plt.xlabel('days')
    plt.legend()
    plt.show()
   
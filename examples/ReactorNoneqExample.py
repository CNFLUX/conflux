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

E_thresh = 1.8
step_thresh = int(E_thresh/binwidth)
runid = 0

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
    
    plt.figure()
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
        plt.plot(days,  rates[key], label=key)
    plt.legend()
    plt.savefig("fission_rate.png")

    spectra_files = {
        "U235": f"{CONFLUX_DB}/example_models/U235_thermal_burst_fission_neutrino.csv",
        "U238": f"{CONFLUX_DB}/example_models/U238_fast_burst_fission_neutrino.csv",
        "Pu239": f"{CONFLUX_DB}/example_models/Pu239_thermal_burst_fission_neutrino.csv",
        "Pu241": f"{CONFLUX_DB}/example_models/Pu241_thermal_burst_fission_neutrino.csv",
    }
    istplist = list(spectra_files.keys())

    spectrum_E_t = {}
    totalflux = {}

    newrates = {} # fission rates that is interpolated from customized fission data
    spectrum_day = {}
    summed_flux = []

    # prepare cumulative fission yield calculation for comparison
    fissile_istp = {}
    fissile_istp["U235"] = FissionIstp(92, 235, 0., DB='JEFF', IFPY=False)
    fissile_istp["U238"] = FissionIstp(92, 238, 0.4, DB='JEFF', IFPY=False)
    fissile_istp["Pu239"] = FissionIstp(94, 239, 0., DB='JEFF', IFPY=False)
    fissile_istp["Pu241"] = FissionIstp(94, 241, 0., DB='JEFF', IFPY=False)
    totalflux_cumu = {}
    spectrum_day_cumu = {}

    # days = np.array([2, 450])
    fig3, ax3 = plt.subplots()
    for j in range(len(istplist)):
        key = istplist[j]

        # loading the S(E,t) data
        df_spectrum = pd.read_csv(spectra_files[key])
        beta_time = df_spectrum.iloc[:, 0].values
        energies = df_spectrum.columns[1:].astype(float).values
        spectrum_E_t[key] = df_spectrum.iloc[:, 1:].values
 
        totalflux[key] = []
        newrates[key] = [] 

        spectrum_day[key] = [] # create a list of spectra over all days

        fissile_istp[key].LoadFissionDB(DB='JEFF')
        fissile_istp[key].LoadCorrelation(DB='JEFF')
        fissile_istp[key].CalcBetaSpectra(betaSpectraDB)
        totalflux_cumu[key] = []
        spectrum_day_cumu[key] = fissile_istp[key].spectrum

        fig1, ax1 = plt.subplots() #  the figure of neutrino flux over time of this figure.
        fig2, ax2 = plt.subplots() #  the figure of neutrino spectrum per fission over time of this figure.
        fig6, ax6 = plt.subplots()

        # rates[key] = np.ones(len(fission_time))
        
        for i in range(len(days)):
            timex = days[i]*3600*24

            # Core function to convolve the neutrino spectrum and the fission rate over a long period 
            spectrum = convolve_via_s_log(timex,                                     
                                        fission_time=fission_time, 
                                        fission_rate=rates[key], 
                                        beta_time=beta_time, 
                                        beta_spec=spectrum_E_t[key], 
                                        N=1000, 
                                        )
            
            spectrum_day[key].append(spectrum) # append the spectrum list of the specified fissile isotope with the spectrum of current day
            totalflux[key].append(sum(spectrum[step_thresh:])*binwidth) # append the neutrino flux list with the flux of current day
            ax1.plot(xbins, spectrum, label = f"{days[i]} days") # plot this day's flux
            # ax2.plot(xbins, spectrum/rates[key], label = f"{days[i]} days")

            newrates[key] = (fission_rate_t(timex, fission_time, rates[key]))
            spectrum_day_cumu[key] = newrates[key]*spectrum_day_cumu[key]
            ax6.plot(xbins, (spectrum-spectrum_day_cumu[key])/spectrum_day_cumu[key], label = f"{days[i]} days")
            totalflux_cumu[key].append(sum(spectrum_day_cumu[key][step_thresh:])*binwidth)
        ax1.plot(xbins, spectrum_day_cumu[key], label = f"cumulative")

        ax1.legend()
        ax1.set_xlabel("energy (MeV)")
        ax1.set_ylabel("neutrinos/MeV")
        fig1.savefig(f"{key}_non_eq_flux_{runid}.png")
        # ax2.legend()
        # ax2.set_xlabel("energy (MeV)")
        # ax2.set_ylabel("neutrinos/MeV/fission")
        # fig2.savefig(f"{key}_non_eq_spectrum_{runid}.png")

        ax3.plot(days, totalflux[key], label = f"{key}")
        # ax3.plot(days, totalflux_cumu[key], label = f"{key}_cumu")
        ax6.legend()
        ax6.set_ylim([-0.3, 0])
        ax6.set_xlabel("energy (MeV)")
        ax6.set_ylabel("residual")
        fig6.savefig(f"{key}_non_eq_spec_residual_{runid}.png")

    ax3.set_xlabel("time (days)")
    ax3.legend()
    fig3.savefig(f"non_eq_summed_flux_{runid}.png")

    totalspect = []
    totalspect_cumu = []
    
    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    for i in range(len(days)):
        spectrum = np.zeros(len(xbins))
        spectrum_cumu = np.zeros(len(xbins))

        for key, item in spectra_files.items():
            spectrum += spectrum_day[key][i]
            spectrum_cumu += spectrum_day_cumu[key][i]
        totalspect.append(sum(spectrum[step_thresh:])*binwidth)
        totalspect_cumu.append(sum(spectrum_cumu[step_thresh:])*binwidth)
        ax5.plot(xbins, spectrum, label=f"{days[i]} days")
        ax5.plot(xbins, spectrum_cumu, label=f"{days[i]} days cumu")
    ax4.plot(days, totalspect, label=f"{days[i]} days")
    ax4.plot(days, totalspect_cumu, label=f"{days[i]} days cumu")
    ax4.set_xlabel('E (MeV)')
    ax4.legend()
    fig4.savefig("total_summed_flux_non_eq.png")

    ax5.set_xlabel('time (days)')
    ax5.legend()
    fig5.savefig("total_flux_non_eq.png")
   
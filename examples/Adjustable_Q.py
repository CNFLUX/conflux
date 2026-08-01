from conflux.FPYEngine import FissionIstp
from conflux.BetaEngine import BetaEngine
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    
    e = np.arange(0, 14., 0.1)
    
    BetaSpec = BetaEngine(xbins = e)
    BetaSpec.LoadBetaDB()
    BetaSpec.CalcBetaSpectra(branchErange=[10, 20])
    
    U235 = FissionIstp(92, 235, Ei = 0, DB = "JEFF")
    U235.LoadFissionDB()
    U235.LoadCorrelation()
    U235.CalcBetaSpectra(BetaSpec)
    
    betaFPY = set(U235.FPYlist.keys()).intersection(BetaSpec.istplist.keys())
    
    for i in betaFPY:
        plt.plot(e, U235.FPYlist[i].y * BetaSpec.istplist[i].spectrum, color = "black")
        plt.fill_between(e, U235.FPYlist[i].y * BetaSpec.istplist[i].spectrum - U235.FPYlist[i].y * BetaSpec.istplist[i].uncertainty,
                         U235.FPYlist[i].y * BetaSpec.istplist[i].spectrum + U235.FPYlist[i].y * BetaSpec.istplist[i].uncertainty, color = "black",
                         alpha = 0.4)
    plt.plot(e, U235.spectrum, color = "red", label = "Total Spectrum")
    plt.fill_between(e, U235.spectrum - U235.uncertainty, U235.spectrum + U235.uncertainty, color = "red",  alpha = 0.4)
    plt.plot([], [], color='black', label='individual branches')
    plt.ylabel(r"$\nu_e$/MeV/fission")
    plt.xlabel("Energy (MeV)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("U235_spectrum.png")
    plt.clf()
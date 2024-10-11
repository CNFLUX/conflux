from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import os
import matplotlib.pyplot as plt
import numpy as np



#One thing that CONFLUX is really good at is making cross database calculations
#As Nuclear data is ever evolving, with better observations, it is good to cross-check
#Various databases for accuracy of data, as well as to check for discrepencies in data.
if __name__ == "__main__":
    #Initially, I will load up a U235 Fission Isotope, one with JEFF, and one with ENDF. I will
    #Then carry out my calculation in the same way. I will define an energy range, load up a 
    #Beta Engine, calculate the BetaSpectra shape of the Fission Products inside the Beta Engine,
    #and then Calculate the Total Beta spectrum for both U235 isotopes.

    U235E = FissionIstp(92, 235, Ei=0, DB='ENDF')
    U235J = FissionIstp(92, 235, Ei=0, DB='JEFF')
    e = np.arange(0,20,0.1)
    BetaSpectraDB = BetaEngine(xbins = e)
    BetaSpectraDB.CalcBetaSpectra(nu_spectrum=False)
    U235E.CalcBetaSpectra(BetaSpectraDB)
    U235J.CalcBetaSpectra(BetaSpectraDB)

    #Notice, I didn't call the LoadFissionDB() method for either U235J or U235E. when initializing
    #a Fission Isotope, LoadFissionDB() is automatically called at the end of the initialization.
    #A nifty little line saver, though be careful that you remember that this is the case.


    #Now that I've calculated the Beta Spectrum for both U235 isotopes, I can go ahead and pull the
    #Spectral information directly from them to plot. Notice, i don't actually need to add my 
    #isotopes into a SummationEngine to pull isotopic information out of individual Fission Isotopes.

    #Here, I calculate the %-difference between ENDF and JEFF
    summedDiff = []
    for i in range(len(e)):
        diff =  (U235E.spectrum[i] - U235J.spectrum[i])/U235E.spectrum[i]
        summedDiff.append(diff)

    #And finally, I plot out both the individual spectrum for each Database, as well as the 
    #%-difference between the calculated spectrum for both DBs
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.fill_between(e, U235E.spectrum + U235E.uncertainty, U235E.spectrum - U235E.uncertainty
                     , color="red", alpha = 0.4, label='ENDF')
    ax1.fill_between(e, U235J.spectrum + U235J.uncertainty, U235J.spectrum - U235J.uncertainty
                     , color="blue", alpha = 0.4, label='JEFF')
    ax2.set(xlabel = "E (in MeV)", ylabel = r"$({\phi}_{ENDF}-{\phi}_{JEFF})/{\phi}_{ENDF}$")
    ax1.set(ylabel = r"${\nu}_e/MeV/fission$")
    ax1.legend()
    ax2.plot(e, summedDiff)
    plt.savefig("JEFFvENDF")

    ax1.set_yscale("log")
    plt.savefig("JEFFvENDF_log")
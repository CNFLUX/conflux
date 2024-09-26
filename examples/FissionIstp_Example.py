from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp, FissionModel
import csv
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":


    # This is an example for the Fission Isotope class inside the FPYEngine.py file
    # The spectrum of the fission isotope must be calculated before adding it to the reactor model.

    #Initially, I set an energy range for my calculation. In this case, 0 - 10 MeV with 100 keV bins.

    e = np.arange(0., 10., 0.1)

    #Here, I've initialized U-235. Note, my atomic number comes before my mass number, and I need to 
    #Specify an incident neutron energy, which is the energy of the neutron that causes the fission
    #Ei can be 0 (thermal), 0.5 (fast), or 14 (relativistic). Note, if I use JEFF, then fast neutrons 
    #Are defined with Ei=0.4 instead of Ei=0.5
    U235 = FissionIstp(92,235, Ei=0)
    
    #After I've initialized my Fission isotope, I can go ahead and load in the correlation and covariance
    #Information of the isotope.
    U235.LoadCorrelation()
    U235.LoadCovariance()

    #One thing to note, I can calculate the spectrum of each beta branch that comes out of this fission
    #isotope. However, in order to do that, I will need a Beta Database, which I can create by intializing
    #A Beta Engine. For more information about that, please look at BetaEngine_Example.py

    #Initialize the beta engine with the energy range I've defined above, and calculate the beta spectra
    #in the given branch Erange.
    BetaEngineDB = BetaEngine(xbins = e)
    BetaEngineDB.CalcBetaSpectra(branchErange=[0,10])

    #Next, I use this BetaEngine to calculate the Spectra for U235. By default, the calculation will
    #Only look at Cumulative fission products. However, one can also calculate the independant fission
    #products in a given time range by passing integer values to the "ifp_begin" and "ifp_end" variables.
    U235.CalcBetaSpectra(BetaEngineDB, processMissing=True)

    #Lastly, I will go ahead and plot each individual beta branch in this fission isotope
    fig = plt.plot()
    for FPZAI in set(U235.FPYlist.keys()).intersection(BetaEngineDB.istplist.keys()):
        individualBranch = U235.FPYlist[FPZAI].y * BetaEngineDB.istplist[FPZAI].spectrum
        plt.plot(e, individualBranch)
    plt.yscale('log')
    plt.savefig("IndividualBranches.png")
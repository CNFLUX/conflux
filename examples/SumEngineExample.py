from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
import numpy as np
if __name__  == "__main__":

    #---------------------------------------------------
    #Full Reactor Spectrum Using U235

    #Set the energy scale of the reaction we want to run. In this case,
    #the calculation is being set from 0 MeV-15 MeV with 100keV bins

    e = np.arange(0, 15., 0.1)

    xbins = e

    # # Calculate beta spectra of all beta unstable isotopes
    betaSpectraDB = BetaEngine(xbins=xbins)
    filename = "beta_spectra.csv"
    try:
        with open(filename, "r") as file:
            betaSpectraDB.LoadFile(filename)
    except FileNotFoundError:
        print("File not found. Creating the file.")
        betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
        betaSpectraDB.SaveToFile(filename)
    print("File created.")
        
    #First, load the BetaEngine and calculate the spectral shapes of all fission branches. This is
    #now our Beta Spectra Database, that we will input into our fission model. 
    # betaSpectraDB = BetaEngine(xbins = e)
    # betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
    
    #initialize the Isotopes that you would like to use in your
    #Reaction and load them into the Database. Here, remember that the
    #first parameter is the Z#, the second is the A#, and the third is
    #Ei, or the incident energy of the neutron that causes fission (0 - Thermal, .5 - Fast, 14 - relativistic)
    U235 = FissionIstp(92,235, Ei= 0)
    U235.LoadFissionDB()
    
    #Load up the Correlation data, and calculate the Covariance matrix (both are done by calling LoadCorrelation).
    #Must be done after loading up the Fission Data
    U235.LoadCorrelation()
    
    #I will quickly calculate the beta spectrum of U235. I have to calculate the U235 spectrum before calculating the
    #Total spectrum. If I add multiple isotopes to the Summation Engine, I will need each isotopes individual
    #Spectrum to be able to calcualte the total spectrum. I will get an error if I do not do this step before
    #Trying to add the isotope to the SummationEngine
    U235.CalcBetaSpectra(betaSpectraDB=betaSpectraDB, processMissing=True)

    #Finally, we are going to load up a summation engine, with the Beta Spectra Database that we just
    #Calculated as an input. As a default, the Summation engine calculates the neutrino specturm,
    #But there is an option to select the beta-spectrum instead.
    SummationEngine = SumEngine(betaSpectraDB=betaSpectraDB)

    #I will also go ahead and load in the Fission Isotope I initialized at the beginning of this example into
    #The Summation Engine. Note the parameters, I add the fission Isotope, the name of the fission Isotope, as
    #Well as how many times I've added that isotope into our summation engine. In this case, I've added it 3 times
    SummationEngine.AddFissionIstp(U235, "U235", count = 1)

    #Finally, I will be calculating the neutrino spectrum for a pure U235 reactor by calling the CalcReactorSpectrum()
    #Function. 
    SummationEngine.CalcReactorSpectrum()

    #And here I will pull out the uncertainties and spectrum from my reactor engine, to plot the spectrum and uncertainties.
    spectrum = SummationEngine.spectrum
    uncertainty = SummationEngine.uncertainty
    FPY_uncertainty = SummationEngine.yieldUnc
    model_uncertainty = SummationEngine.modelUnc

    #Draw the SummationEngineing spectra
    fig = plt.figure()

    #Here, I am plotting all the individual beta branches that make up the totality of the
    #Reactor spectrum. It is a good visualization of what the "ab-initio" method does when
    #Making it's neutrino prediction
    loopcount = 0
    for FPZAI in set(SummationEngine.FPYlist.keys()).intersection(betaSpectraDB.istplist.keys()):
        individualBranch = SummationEngine.FPYlist[FPZAI].y * betaSpectraDB.istplist[FPZAI].spectrum
        label = "Individual FP spectra" if loopcount == 0 else "_nolegend_"
        plt.plot(e, individualBranch, linestyle="--", color="grey", label=label)
        loopcount += 1


    #This is the plotting of the total spectrum
    plt.errorbar(e, SummationEngine.spectrum, yerr=uncertainty, label="Summed neurtrino spectrum")
    plt.ylim([1e-4, 10])
    plt.yscale('log')
    plt.legend()
    plt.xlabel("E (MeV)")
    plt.ylabel("neutrinos/MeV/Fission")
    fig.savefig("summation_neutrino_spectrum.pdf")

    #Draw the relative uncertainty 
    rel_err_tot = uncertainty/spectrum*100
    rel_err_mod = model_uncertainty/spectrum*100
    rel_err_fpy = FPY_uncertainty/spectrum*100
    fig = plt.figure()
    # plt.plot(e, rel_err_tot, label='total uncertainty')
    # plt.plot(e, rel_err_mod, label='beta uncertainty')
    # plt.plot(e, rel_err_fpy, label='FPY uncertainty')
    plt.fill_between(e, rel_err_tot, -rel_err_tot, label='total uncertainty', alpha=0.4)
    plt.fill_between(e, rel_err_mod, -rel_err_mod, label='beta model uncertainty', alpha=0.4)
    plt.fill_between(e, rel_err_fpy, -rel_err_fpy, label='FPY uncertainty', alpha=0.4)
    plt.ylim([-100, 100])
    plt.legend()
    plt.xlabel("E (MeV)")
    plt.ylabel("Relative error (%)")
    fig.savefig("relative_summation_uncertainty.pdf")

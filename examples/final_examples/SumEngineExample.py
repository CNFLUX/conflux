from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
import numpy as np
if __name__  == "__main__":

    #---------------------------------------------------
    #Full Reactor Spectrum Using U235

    #Set the energy scale of the reaction we want to run. In this case,
    #the calculation is being set from 0 MeV-15 MeV with 100keV bins

    e = np.arange(0, 15., 0.1)

    #initialize the Isotopes that you would like to use in your
    #Reaction and load them into the Database. Here, remember that the
    #first parameter is the Z#, the second is the A#, and the third is
    #Ei, or the incident energy of the neutron that causes fission (0 - Thermal, .5 - Fast, 14 - relativistic)
    U235 = FissionIstp(92,235, Ei= 0)
    U235.LoadFissionDB()
    #Load up the Correlation data, and calculate the Covariance matrix (both are done by calling LoadCorrelation).
    #Must be done after loading up the Fission Data
    U235.LoadCorrelation()

    #Next, load in a BetaEngine and calculate the spectral shapes of all fission branches. This is
    #now our Beta Spectra Database, that we will input into our fission model. Also note, branchErange
    #is the range of branch endpoint energies we want to look at. For this example, 
    #We are setting the endpoint energy range from 0MeV to 15 MeV
    betaSpectraDB = BetaEngine(xbins = e)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 15.0])

    #I will quickly calculate the beta spectrum of U235. I have to calculate the U235 spectrum before calculating the
    #Total spectrum. If I add multiple isotopes to the Summation Engine, I will need each isotopes individual
    #Spectrum to be able to calcualte the total spectrum. I will get an error if I do not do this step before
    #Trying to add the isotope to the summation_engine
    U235.CalcBetaSpectra(betaSpectraDB=betaSpectraDB)

    #Finally, we are going to load up a summation engine, with the Beta Spectra Database that we just
    #Calculated as an input. As a default, the Summation engine calculates the neutrino specturm,
    #But there is an option to select the beta-spectrum instead.
    summation_engine = SumEngine(betaSpectraDB= betaSpectraDB, nu_spectrum=True)

    #I will also go ahead and load in the Fission Isotope I initialized at the beginning of this example into
    #The Summation Engine. Note the parameters, I add the fission Isotope, the name of the fission Isotope, as
    #Well as how many times I've added that isotope into our summation engine. In this case, I've added it 3 times
    summation_engine.AddFissionIstp(U235, "U235", count = 1)

    #Finally, I will be calculating the neutrino spectrum for a pure U235 reactor by calling the CalcReactorSpectrum()
    #Function. 
    summation_engine.CalcReactorSpectrum()

    #And here I will pull out the uncertainties and spectrum from my reactor engine, to plot the spectrum and uncertainties.
    spectrum = summation_engine.spectrum
    uncertainty = summation_engine.uncertainty

    #Draw the summation_engineing spectra
    fig = plt.figure()


    #Here, I am plotting all the individual beta branches that make up the totality of the
    #Reactor spectrum. It is a good visualization of what the "ab-initio" method does when
    #Making it's neutrino prediction
    for FPZAI in set(summation_engine.FPYlist.keys()).intersection(betaSpectraDB.istplist.keys()):
        individualBranch = summation_engine.FPYlist[FPZAI].y * betaSpectraDB.istplist[FPZAI].spectrum
        plt.plot(e, individualBranch)



    #This is the plotting of the spectra that made up the total example summation spectrum
    plt.errorbar(e, summation_engine.spectrum, yerr=uncertainty)
    plt.plot(e, spectrum, 'b--')
    plt.ylim([1e-8, 10])
    plt.yscale('log')
    plt.xlabel("E (MeV)", fontsize= 18)
    plt.ylabel("neutrinos/MeV/Fission", fontsize = 18)
    # plt.title("U235 spectrum and summed fission product neutrino spectra", fontsize=20)
    fig.savefig("U235_Spectra_calculated.png")
    
    #This is the plotting of the total spectrum
    fig = plt.figure()
    plt.errorbar(e, summation_engine.spectrum, yerr=uncertainty)
    plt.plot(e, spectrum, 'b--')
    plt.xlabel("E (MeV)", fontsize= 18)
    plt.ylabel("neutrinos/MeV/Fission", fontsize = 18)
    plt.title("U235 Neutrino Spectrum", fontsize=20)
    fig.savefig("U235_Neutrino_Spectrum.png")


    # This is to plot the contribution of each type of uncertainties
    total_err = summation_engine.uncertainty
    model_err = summation_engine.modelUnc
    yield_err = summation_engine.yieldUnc
    
    fig, ax = plt.subplots()
    ax.set_xlim([0, 10])
    ax.set_ylim([-30, 30])
    ax.set(xlabel='E (MeV)', ylabel='relative error (%)')
    ax.fill_between(e, total_err/spectrum*100, -total_err/spectrum*100, label="total uncertainty", alpha=.5)
    ax.fill_between(e, yield_err/spectrum*100, -yield_err/spectrum*100, label="fission product uncertainty", alpha=.5)
    ax.fill_between(e, model_err/spectrum*100, -model_err/spectrum*100, label="beta model uncertainty", alpha=.5)
    
    ax.legend()
    fig.savefig("U235_spectrum_uncertainty.png")
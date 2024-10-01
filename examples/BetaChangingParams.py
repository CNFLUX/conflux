from conflux.BetaEngine import BetaIstp, BetaEngine
from conflux.bsg.SpectralFunctions import fermi_function, phase_space
from conflux.bsg.Constants import ELECTRON_MASS_MEV
import numpy as np
import matplotlib.pyplot as plt
import copy

#here, I show an interesting feature of CONFLUX; editing various theretical parameters of
#a beta isotope


#I've added a custom function to calculate a beta spectrum. I will be using this later in the
#program to load in a custom function to calculate the beta spectrum instead of the default 
#function with all the shape corrections.

def ferm(ebeta, p, numass = 0):
    result = 0.
    W = ebeta/ELECTRON_MASS_MEV + 1
    result = fermi_function(W, **p) * phase_space(W, **p, numass=numass)

    return result

if __name__ == "__main__":
    #First, initialize energy range
    e = np.arange(0., 0.5, .001)

    #Next, I am going to initialize a beta engine with one Beta Isotope, and
    #Pull the Beta Isotope from the engine after I've loaded the beta data. 
    #This is to ensure that the BetaIstp object that I work with has all isotopic
    #Information from the ENSDFbetaDB2

    #One thing I want to point out, even though I've only added one isotope to the 
    #BetaEngine, when I run CalcBetaSpectrum, it tells me that two isotopes are being
    #Calculated. This is because Cs141 has 2 isomeric states to calculate. We differentiate
    #The two with the last digit of the isotopes FPZAI number. Hence, the original isotope
    #is 551410, and the first isomeric state is 551411
    BetaDB = BetaEngine(inputlist = [551410], xbins = e)
    BetaDB.LoadBetaDB()
    BetaDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0,8])
    Cs141 = BetaDB.istplist[551410]

    #I also go ahead and save the original spectrum and uncertainty before making any edits
    #To my theoretical model.
    origSpec = (Cs141.spectrum)
    origUnc = (Cs141.uncertainty)
    
    #Next, I create a new branch with a slightly higher endpoint energy than the 
    #highest observed branch, with some given parameters

    #In this case, 'fraction' is the fractional contribution this branch has on the total Cs141 
    # spectrum, 'sigma_frac' is the uncertainty in the fractional contribution, and 'sigma_E0'
    # is the uncertainty in the endpoint energy of this branch.
    parameters = {'fraction' : 0.1, 'sigma_E0' : 0.01, 'sigma_frac' : 0.01}
    Cs141.EditBranch(defaultE0 = 5.244, kwargs =  parameters)
    #I also go ahead and recalculate the covariances of all branches, since I just edited the spectrum
    #To add a branch that did not previously exist.
    Cs141.CalcCovariance()

    #And here I recalculate the total spectrum with the extra branch that I just defined. I also save the
    #resulting spectrum and uncertainties.
    Cs141.SumSpectra(nu_spectrum = True, branchErange=[0., 8])
    EndSpec = Cs141.spectrum
    EndUnc = Cs141.uncertainty
    
    # #So, I just added a branch with a different endpoint energy, let me do something else. I can also edit
    # #The forbiddenness of a branch. First, I will remove the branch that I added above, then iterate through
    # #All branches and change their forbiddenness to 1. After, since I've edited a few branches and removed one,
    # #I will recalculate the covariances for my spectrum and then save the resulting spectrum
    
    #Here I remove the branch I just added
    Cs141.branches.pop(5.244)
    #And here I change All the branches forbiddenness to 0. 
    for i in Cs141.branches:
        if Cs141.branches[i].forbiddenness != 1:
            Cs141.EditBranch(defaultE0 = i, forbiddenness = 0)
    Cs141.CalcCovariance()
    Cs141.SumSpectra(nu_spectrum=True, branchErange=[0,8])

    forbidSpec = Cs141.spectrum
    forbidUnc = Cs141.uncertainty

    #finally, I'm going to recalculate the beta spectrum without any spectral shape corrections
    #Only calculating it with the fermi function

    #First, I will go ahead and pull all the isotopic information from my ENSDFbetaDB2 again, and
    # create a fresh copy of the Cs141 isotope from the DB, since I want to start this last calculation
    # off with a fresh BetaIstp That has not been edited.
    BetaDB.LoadBetaDB()
    BetaDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0,8])
    Cs141 = BetaDB.istplist[551410]

    for i in Cs141.branches:
        Cs141.EditBranch(defaultE0 = i, custom_func = ferm)
    #Then I go ahead and Calculate my beta spectrum, since I have edited branches, and save the resulting
    #Spectrum and uncertainties
    Cs141.SumSpectra(nu_spectrum=True, branchErange=[0,8])
    fermiSpec = (Cs141.spectrum)
    fermiUnc = (Cs141.uncertainty)

    #And finally, I can go ahead and plot all of this, and save it in a file

    fig = plt.figure()

    plt.errorbar(e, EndSpec, yerr=EndUnc, label="Endpoint edited Spectrum")
    plt.errorbar(e, forbidSpec, yerr=forbidUnc, label="Forbiddenness edited")
    plt.errorbar(e, fermiSpec, yerr = fermiUnc, label="beta function edited")
    plt.errorbar(e, origSpec,yerr=origUnc, label="Vanilla Spectrum")
    plt.ylabel("neutrino Spectrum")
    plt.xlabel("Energy (in MeV)")
    plt.legend()
    plt.savefig("Beta_spectrum_changing_params.png")
    plt.yscale("log")
    plt.savefig("Beta_Spectrum_changing_params_log.png")
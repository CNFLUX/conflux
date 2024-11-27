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
    
    result = np.nan_to_num(result, nan=0.0)

    return result

def safe_division(a, b):
    result = np.where(b != 0, a / b, 0)
    return result

if __name__ == "__main__":
    #First, initialize energy range
    e = np.arange(0., 8, 0.1)

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
    origSpec =  BetaDB.istplist[551410].spectrum
    origUnc =  BetaDB.istplist[551410].uncertainty
    
    #Next, I create a new branch with a slightly higher endpoint energy than the 
    #highest observed branch, with some given parameters

    #In this case, 'fraction' is the fractional contribution this branch has on the total Cs141 
    # spectrum, 'sigma_frac' is the uncertainty in the fractional contribution, and 'sigma_E0'
    # is the uncertainty in the endpoint energy of this branch.
    parameters = {'fraction' : 0.1, 'sigma_E0' : 0.01, 'sigma_frac' : 0.01}
    parameters1 = {'E0' : 6.2}
    # When I gave an E0 argument to a defaultE0, that E0 will replace the defaultE0
    Cs141.EditBranch(defaultE0 = 5.206, E0=6.206)
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
    
    
    # finally, I'm going to recalculate the beta spectrum without any spectral shape corrections
    # Only calculating it with the fermi function

    # First, I will go ahead and undo all of the changes to the forbiddenness values for the 
    # Beta Branches, and then edit the branches so that they take in a custom beta spectra
    # function as I've defined above (The fermi function mutliplied by the phase space, so without
    # any theoretical corrections)
    BetaDB3 = BetaEngine(inputlist = [551410], xbins = e)
    BetaDB3.LoadBetaDB()
    Cs141_3 =BetaDB3.istplist[551410]
    
    # switch all allowed transitions to forbidden transitions
    for i in Cs141_3.branches:
        if Cs141_3.branches[i].forbiddenness == 0:
            Cs141_3.EditBranch(defaultE0 = i, forbiddenness = 1) 
        
    # recalculate the beta spectra 
    BetaDB3.CalcBetaSpectra(nu_spectrum=True)    
    
    forbidSpec = Cs141_3.spectrum
    forbidUnc = Cs141_3.uncertainty

    # Switch the default beta shape function with a custom function I've defined up top.
    # First, I will go ahead and pull all the isotopic information from my ENSDFbetaDB2 again, and
    # create a fresh copy of the Cs141 isotope from the DB, since I want to start this last calculation
    # off with a fresh BetaIstp That has not been edited.
    BetaDB2 = BetaEngine(inputlist = [551410], xbins = e, custom_func=ferm)
    BetaDB2.LoadBetaDB()
    BetaDB2.CalcBetaSpectra(nu_spectrum=True, branchErange=[0,8])
    Cs141_2 = BetaDB2.istplist[551410]

    # for i in Cs141.branches:
    #     Cs141_2.EditBranch(defaultE0 = i, custom_func = ferm)
    # Cs141_2.Display()
    # #Then I go ahead and Calculate my beta spectrum, since I have edited branches, and save the resulting
    # #Spectrum and uncertainties
    # Cs141_2.SumSpectra(nu_spectrum=False, branchErange=[0,8])
    fermiSpec = (Cs141_2.spectrum)
    fermiUnc = (Cs141_2.uncertainty)

    #And finally, I can go ahead and plot all of this, and save it in a file

    fig = plt.figure()
    plt.plot(e, safe_division(forbidSpec-origSpec, origSpec)*100, label="Forbiddenness edited")
    plt.ylabel("residual (%)")
    plt.xlabel("Energy (MeV)")
    plt.legend()
    plt.savefig("Beta_spectrum_changing_forbid.pdf")
    
    fig = plt.figure()
    plt.plot(e, (EndSpec), label="Endpoint edited")
    plt.plot(e, (fermiSpec), label="Customized fermi function")
    plt.plot(e, (origSpec), label="CONFLUX default output")
    plt.yscale("log")
    plt.ylim([1e-4, 1])
    plt.ylabel("neutrino/MeV")
    plt.xlabel("Energy (MeV)")
    plt.legend()
    plt.savefig("Beta_spectrum_changing_params.pdf")


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:11:30 2024

@author: Xianyi Zhang (LLNL)
"""

from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
import numpy as np

if __name__  == "__main__":
    # setup the spectrum energy range and binning
    e = np.arange(0, 15., 0.1)
    
    # Calculate the neutrino spectrum by setting the missing branches to be 3 
    # for each missing isotopes. The branches are set to be rate = 1/3 fraction 
    # each with a 1/3*Q, 2/3*Q, and Q value as the end point energy, by setting
    # `missingBranch = 3`, which is also the default value.
    betaSpectraDB = BetaEngine(xbins = e)
    betaSpectraDB.LoadBetaDB(missingBranch=3)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
    
    # calculate the specrtra with an alternative example where only one branch
    # is set for missing isotopes
    betaSpectraDB1 = BetaEngine(xbins = e)
    betaSpectraDB1.LoadBetaDB(missingBranch=1)
    betaSpectraDB1.CalcBetaSpectra(nu_spectrum=True)
    
    # Load the fission products of 235U with thermal neutron fission "Ei = 0"
    U235 = FissionIstp(92,235, Ei= 0)
    U235.LoadFissionDB()
    # Load the default correlation matrix saved in the $CONFLUX_DB, it will 
    # calculate the covariance matrix automatically for you.
    U235.LoadCorrelation()
    # Calculate the neutrino spectrum by calculating spectrum 
    U235.CalcBetaSpectra(betaSpectraDB=betaSpectraDB, processMissing=True)
    # Summing all fissile isotope spectra. In this case, only one fission of 
    # 235U.
    # Picking the 3-branch missing addition in the neutrino spectrum calculation
    SummationEngine = SumEngine(betaSpectraDB=betaSpectraDB)
    SummationEngine.AddFissionIstp(U235, "U235", count = 1)
    SummationEngine.CalcReactorSpectrum()
    #And here I will pull out the uncertainties and spectrum from my reactor engine, to plot the spectrum and uncertainties.
    spectrum = SummationEngine.spectrum
    uncertainty = SummationEngine.uncertainty
    
    # U235 spectrum with missing info added with 1-branch assumption
    U2351 = FissionIstp(92,235, Ei= 0)
    U2351.LoadFissionDB()
    U2351.LoadCorrelation()
    U2351.CalcBetaSpectra(betaSpectraDB=betaSpectraDB1, processMissing=True)
    SummationEngine1 = SumEngine(betaSpectraDB=betaSpectraDB1)
    SummationEngine1.AddFissionIstp(U2351, "U235_1", count = 1)
    SummationEngine1.CalcReactorSpectrum()

    #And here I will pull out the uncertainties and spectrum from my reactor engine, to plot the spectrum and uncertainties.
    spectrum1 = SummationEngine1.spectrum
    uncertainty1 = SummationEngine1.uncertainty
  
    # U235 spectrum with no missing info added
    # in this case, picking beta spectrum db as betaSpectraDB1 or betaSpectraDB doesn't matter anymore
    U2352 = FissionIstp(92,235, Ei= 0)
    U2352.LoadFissionDB()
    U2352.LoadCorrelation()
    U2352.CalcBetaSpectra(betaSpectraDB=betaSpectraDB1, processMissing=False)
    SummationEngine1 = SumEngine(betaSpectraDB=betaSpectraDB1)
    SummationEngine1.AddFissionIstp(U2352, "U235_2", count = 1)
    SummationEngine1.CalcReactorSpectrum()
    spectrum2 = SummationEngine1.spectrum
    uncertainty2 = SummationEngine1.uncertainty
    
    #Draw the SummationEngineing spectra
    fig = plt.figure()

    #This is the plotting of the total spectrum
    plt.plot(e, spectrum2, label='spectrum w/o missing info')
    plt.errorbar(e, spectrum1-spectrum2, label='1-branch missing contribution')
    plt.errorbar(e, spectrum-spectrum2, label='3-branch missing contribution')
    # plt.xlim([1.8, 10])
    plt.ylim([1e-4, 10])
    plt.yscale('log')
    plt.xlabel("E (MeV)")
    plt.legend()
    plt.ylabel("relative contribution of missing info")
    fig.savefig("summation_neutrino_spectrum.pdf")


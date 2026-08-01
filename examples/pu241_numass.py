#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 20:22:13 2024

@author: zhang39
"""
import numpy as np
import matplotlib.pyplot as plt

from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp, FissionModel


if __name__ == "__main__":
    x = np.arange(0, 0.04, 0.0001)

    testlist = [942410]
    testEngine = BetaEngine(xbins=x, inputlist=testlist, numass=8e-3, mixing=0.5)
    testEngine.istplist[942410].EditBranch(0.021, forbiddenness = 1)
    testEngine.CalcBetaSpectra(nu_spectrum=False, GSF=False)

    y = testEngine.istplist[942410].spectrum
    yerr = testEngine.istplist[942410].uncertainty
    
    testEngine.istplist[942410].Display()

    testEngine0 = BetaEngine(xbins=x, inputlist=testlist)
    testEngine0.istplist[942410].EditBranch(0.021, forbiddenness = 1)
    testEngine0.CalcBetaSpectra(nu_spectrum=False, GSF=False)


    y0 = testEngine0.istplist[942410].spectrum
    y0err = testEngine0.istplist[942410].uncertainty

    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.errorbar(x/1e-3, y/1e3, label="8 keV sterile neutrino", yerr=yerr)
    ax1.errorbar(x/1e-3, y0/1e3, label="zero neutrino mass", yerr=y0err)
    ax1.set(ylabel=("beta/keV") )
    ax1.legend()
    ax2.errorbar(x/1e-3, (y-y0)/1e3, label="8 keV sterile neutrino", yerr=yerr)
    ax2.set(xlabel="E (keV)", ylabel=("difference") )
    ax1.legend()
    plt.savefig("pu241_sterile_neutrino2.pdf")
    
    # U235 = FissionIstp(92,235, Ei=0)
    # U235.LoadCorrelation()
    # U235.LoadCovariance()
    # U235.CalcBetaSpectra(testEngine, processMissing=True)
    # spect1 = U235.spectrum
    # unc1 = U235.uncertainty
    
    # U235.CalcBetaSpectra(testEngine0, processMissing=True)
    # spect0 = U235.spectrum
    # unc0 = U235.uncertainty
    
    # plt.xlabel("E (MeV)")
    # plt.ylabel("beta/MeV")
    # plt.ylim([-0.05, 0.05])
    # plt.errorbar(x, (spect1-spect0)/spect0, label="8 keV sterile neutrino", yerr=0)
    # # plt.errorbar(x, spect0, label="zero neutrino mass", yerr=unc0)
    
    # plt.legend()
    # plt.savefig("Fission_steril_neutrion.pdf")
    
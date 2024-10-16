#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 20:22:13 2024

@author: zhang39
"""
import numpy as np
import matplotlib.pyplot as plt

from conflux.BetaEngine import BetaEngine

if __name__ == "__main__":
    x = np.arange(0, 50e-3, 1e-4)
    binwidth = 1

    testlist = [942410]
    testEngine = BetaEngine(testlist, xbins=x, numass=8e-3, mixing=0.5)
    testEngine.CalcBetaSpectra(nu_spectrum=False, branchErange=[0.0, 20], GSF=False)

    y = testEngine.istplist[942410].spectrum
    yerr = testEngine.istplist[942410].uncertainty
    
    testEngine0 = BetaEngine(testlist, xbins=x)
    testEngine0.CalcBetaSpectra(nu_spectrum=False, branchErange=[0.0, 20], GSF=False)

    y0 = testEngine0.istplist[942410].spectrum
    y0err = testEngine0.istplist[942410].uncertainty
    
    plt.xlabel("E (keV)")
    plt.ylabel("beta/keV")
    plt.errorbar(x/1e-3, y/1e3, label="8 keV sterile neutrino", yerr=yerr)
    plt.errorbar(x/1e-3, y0/1e3, label="zero neutrino mass", yerr=y0err)
    
    
    plt.legend()
    plt.savefig("pu241_sterile_neutrino.pdf")
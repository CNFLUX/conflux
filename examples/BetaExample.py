# An example of the beta engine
# TODO: Descriptions goes here

import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaEngine import BetaEngine

if __name__ == "__main__":
    x = np.arange(1.5, 10, 0.05)
    binwidth = 1

    #Initialize a list of beta isotopes you would like to carry the beta calculation for.
    #Note, that all isotopes are ID'd in ZAI format (First two digits are Z[Atomic Number], next 3 are A [Neutron number], and last digit is I [Isomer number])
    testlist = [390960, 390961, 521331, 531371, 922390, 932390]
    #Initialize the beta Engine with the user defined list of isotopes
    testEngine = BetaEngine(testlist, xbins=x)
    testEngine.CalcBetaSpectra(nu_spectrum=False)
    
    #print(testEngine.spectralist[521340])
    index = 521331
    fig = plt.figure()
    #plot out a specific isotope based on it's ZAI number.
    y1 = testEngine.istplist[index].spectrum*binwidth
    plt.errorbar(x, y1, label=str(index)+"_beta") #, yerr=testEngine.istplist[index].uncertainty)
    # print(max(testEngine.istplist[531340].spectrum), sum(testEngine.istplist[531340].spectrum))
    # testEngine.CalcBetaSpectra(nu_spectrum=False)
    # plt.errorbar(x, testEngine.istplist[531340].spectrum, yerr=testEngine.istplist[531340].uncertainty)
    # print(max(testEngine.istplist[531340].spectrum), sum(testEngine.istplist[531340].spectrum))

    #initialize a seconnd beta engine so that you can do both Beta and Neutrino calculations.
    testEngine2 = BetaEngine(testlist, xbins=x)
    testEngine2.CalcBetaSpectra(nu_spectrum=True)
    y2 = testEngine2.istplist[index].spectrum*binwidth
    print("beta/nu", sum(y1))
    plt.xlabel("E (MeV)")
    plt.plot(x, y2, "--", label=str(index)+"_nu") #, yerr=testEngine.istplist[index].uncertainty)
    plt.legend()
    plt.show()
    
    # print(testEngine2.istplist[index].branches)

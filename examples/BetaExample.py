# An example of the beta engine
# TODO: Descriptions goes here

import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaEngine import BetaEngine

if __name__ == "__main__":
    x = np.linspace(0, 20, 200)

    testlist = [390960, 390961, 521340, 531340]
    testEngine = BetaEngine(testlist, xbins=x)
    testEngine.CalcBetaSpectra(nu_spectrum=True)
    
    #print(testEngine.spectralist[521340])

    fig = plt.figure()
    plt.errorbar(x, testEngine.istplist[531340].spectrum, yerr=testEngine.istplist[531340].uncertainty)
    print(max(testEngine.istplist[531340].spectrum), sum(testEngine.istplist[531340].spectrum))
    testEngine.CalcBetaSpectra(nu_spectrum=False)
    plt.errorbar(x, testEngine.istplist[531340].spectrum, yerr=testEngine.istplist[531340].uncertainty)
    print(max(testEngine.istplist[531340].spectrum), sum(testEngine.istplist[531340].spectrum))

    plt.show()
    # fig.savefig("errorbartest.png")

    #testbeta = BetaBranch(1, 3, 1.0, 0, )

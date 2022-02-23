# An example of the beta engine
# TODO: Descriptions goes here

import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaEngine import BetaEngine

if __name__ == "__main__":
    testlist = [390960, 390961, 521340, 531340]
    testEngine = BetaEngine(testlist)
    testEngine.CalcBetaSpectra(nu_spectrum=True)
    #print(testEngine.spectralist[521340])

    fig = plt.figure()
    x = np.linspace(0, 20, 200)
    plt.errorbar(x, testEngine.spectralist[390960], yerr=testEngine.uncertaintylist[390961])
    #plt.draw()
    fig.savefig("errorbartest.png")

    #testbeta = BetaBranch(1, 3, 1.0, 0, )

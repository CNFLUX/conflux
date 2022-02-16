# An example of the beta engine
# TODO: Descriptions goes here

import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaEngine import BetaEngineif __name__ == "__main__":
    testlist = [521340, 531340, 350900]
    testEngine = BetaEngine(testlist)
    testEngine.CalcBetaSpectra(nu_spectrum=True)
    #print(testEngine.spectralist[521340])
    print(testEngine.spectralist[350900])
    print(testEngine.uncertaintylist[350900])

    fig = plt.figure()
    x = np.linspace(0, 20, 200)
    plt.errorbar(x, testEngine.spectralist[350900], yerr=testEngine.uncertaintylist[350900])
    #plt.draw()
    fig.savefig("errorbartest.png")

    #testbeta = BetaBranch(1, 3, 1.0, 0, )

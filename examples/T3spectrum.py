# A tritium beta spectrum generator

import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaEngine import BetaEngine


if __name__ == "__main__":
    testlist = [942410]
    testEngine = BetaEngine(testlist)
    testEngine.CalcBetaSpectra(nu_spectrum=False, binwidths=0.0001, lower=-1.0, thresh=0.0, erange = 0.02)
    #print(testEngine.spectralist[10030])
    x = np.linspace(0., 0.02, 200)
    y = np.array(testEngine.spectralist[942410]/10000)
    print(y.sum())
    fig = plt.figure()
    plt.errorbar(x, y)
    fig.savefig('test/Pu241beta.png')

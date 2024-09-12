# An example of the beta engine
# TODO: Descriptions goes here

import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaEngine import BetaEngine

if __name__ == "__main__":
    x = np.arange(0, 1, 0.005)
    binwidth = 0.005

    testEngine = BetaEngine(xbins=x)
    testEngine.CalcBetaSpectra(nu_spectrum=False, branchErange=[0.0, 0.1])

    istp = 922370 #U237
    testlist = [istp]
    fig = plt.figure()

    testEngine2 = BetaEngine(testlist, xbins=x)
    testEngine2.CalcBetaSpectra(nu_spectrum=False)
    y2 = testEngine2.istplist[istp].spectrum*binwidth
    y2err = testEngine2.istplist[istp].uncertainty*binwidth
    print("beta/nu", sum(y2))
    plt.xlabel("E (keV)")
    plt.errorbar(x/1e-3, y2, label=str(istp)+"_nu", yerr=y2err)
    plt.legend()
    plt.show()

    # print(testEngine2.istplist[index].branches)

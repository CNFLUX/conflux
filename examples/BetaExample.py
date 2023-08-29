# An example of the beta engine
# TODO: Descriptions goes here

import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaEngine import BetaEngine

if __name__ == "__main__":
    x = np.linspace(0, 4, 200)
    binwidth =  1

    testlist = [390960, 390961, 521340, 531340, 922390, 932390]
    testEngine = BetaEngine(testlist, xbins=x)
    testEngine.CalcBetaSpectra(nu_spectrum=False)
    
    #print(testEngine.spectralist[521340])
    index = 922390
    fig = plt.figure()
    y1 = testEngine.istplist[index].spectrum*binwidth
    plt.errorbar(x, y1, label=str(index)+"_beta") #, yerr=testEngine.istplist[index].uncertainty)
    # print(max(testEngine.istplist[531340].spectrum), sum(testEngine.istplist[531340].spectrum))
    # testEngine.CalcBetaSpectra(nu_spectrum=False)
    # plt.errorbar(x, testEngine.istplist[531340].spectrum, yerr=testEngine.istplist[531340].uncertainty)
    # print(max(testEngine.istplist[531340].spectrum), sum(testEngine.istplist[531340].spectrum))

    
    testEngine2 = BetaEngine(testlist, xbins=x)
    testEngine2.CalcBetaSpectra(nu_spectrum=True)
    y2 = testEngine2.istplist[index].spectrum*binwidth
    print("beta/nu", sum(y1)/sum(y2))
    plt.xlabel("E (MeV)")
    plt.plot(x, y2, "--", label=str(index)+"_nu") #, yerr=testEngine.istplist[index].uncertainty)
    plt.legend()
    plt.show()
    
    # print(testEngine2.istplist[index].branches)

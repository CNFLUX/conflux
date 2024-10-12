from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
import numpy as np
if __name__  == "__main__":

# Y96 ,Rb92 ,Cs142 ,Y97 ,Rb93 ,Nb100 ,Cs140 ,Sr95 
    isolist = [390960, 370920, 551420, 390970, 370930, 411001, 551400, 380950]

    e = np.arange(0., 10., 0.1)
    BetaSpecificEngine = BetaEngine(inputlist=isolist, xbins = e)
    BetaSpecificEngine.CalcBetaSpectra()

    fig = plt.figure()

    for i in BetaSpecificEngine.istplist:
        plt.errorbar(e, BetaSpecificEngine.istplist[i].spectrum, 
                     yerr = BetaSpecificEngine.istplist[i].uncertainty, label = BetaSpecificEngine.istplist[i].name)
        
    plt.xlabel("Energy (MeV)")
    plt.ylabel(r"$e^-/MeV/Fission$")
    plt.legend()
    plt.savefig("Specific_Neutrino_Spectrum.pdf")

    #Will need to create individual beta isotopes with the neutrino param set to True
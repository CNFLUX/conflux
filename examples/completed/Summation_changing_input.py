from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Initialize U235 Isotope
    U235 = FissionIstp(92,235)
    U235.LoadFissionDB()
    #Initialize the Fission Model, add the U235 Isotope to it
    model = FissionModel()
    model.AddContribution(isotope=U235, Ei=0, fraction=1)
    #Initialize the summation engine, add the Fission model to it.
    result = SumEngine()
    result.AddModel(model)
    #Run the Summation calculation to calculate the reactor spectrum with
    #Default forbidenness
    betaSpectraDB = BetaEngine(result.FPYlist.keys())
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0,20.0])

    result.CalcReactorSpectrum(betaSpectraDB)

    fig = plt.figure()

    #Energy values from the Spectrum
    sumX = result.xbins

    temp = [330840, 330860, 350900]
    plt.plot(sumX, result.spectrum, label = "total original spectrum")

    #Plot out the three isotopes in temp
    for xi in temp:
       print(xi, result.FPYlist[xi].y)
       plt.plot(sumX, result.betaSpectraList[xi] * result.FPYlist[xi].y, label = str(xi))

    #Initialize the second Summation engine that will
    #Take in our changing input
    otherResult = SumEngine()
    otherResult.AddModel(model)

    betalist = otherResult.FPYlist.keys()
    betaDB = BetaEngine(betalist)


    #Here I take a look at each isotope inside the betalist
    for i in betalist:
        #Check if the isotope is inside the isotope list.
        if i in betaDB.istplist:
            #If the isotope is inside the list, iterate through
            #Each of it's branches and change the forbiddeness to -2
            for j in betaDB.istplist[i].branches:
                frac = betaDB.istplist[i].branches[j].frac
                energy = betaDB.istplist[i].branches[j].E0
                betaDB.istplist[i].EditBranch(E0=energy, fraction=frac, forbiddeness=-2)


    #Run the reactor spectrum calculation on the second summation engine
    betaDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])
    otherResult.CalcReactorSpectrum(betaDB)

    #Plot out the spectrum, as well as branches of note
    plt.plot(sumX, otherResult.spectrum, label = 'total changed spectrum', linestyle="dashed")
    for xi in temp:
        x = str(xi) + " changed"
        plt.plot(sumX, otherResult.betaSpectraList[xi] * otherResult.FPYlist[xi].y, label = x, linestyle="dashed")
    #plt.plot(sumX, otherResult.betaSpectraList[330840] * otherResult.FPYlist[330840].y, label = "330840 - Changed")


    plt.legend()
    #plt.title("Change in Forbiddeness to -10 for 3 isotopes of Arsenic compared to the original")
    plt.ylabel("neutrinos/MeV/Fission")
    plt.xlabel("E (in MeV)")
    fig.savefig("U-235Spec_reg_v4")
    plt.yscale('log')
    fig.savefig("U-235_log_v4")

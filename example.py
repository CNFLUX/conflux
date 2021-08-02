from BetaEngine import BetaEngine
from FPYEngine import FissionModel, FissionIstp
from SumEngine import SumEngine

if __name__  == "__main__":

#---------------------------------------------------
    #Full Reactor Spectrum Using U235


    #initialize the Isotopes that you would like to use in your
    #Reaction and load them into the Database

    U235 = FissionIstp(92,235)
    U235.LoadDB()
    U238 = FissionIstp(92, 238)
    U238.LoadDB()
    Pu240 = FissionIstp(94, 240)
    Pu240.LoadDB()
    Pu241 = FissionIstp(94, 241)
    Pu241.LoadDB()

#Initialize the model that you would like to work with.

    modelFullSpec = FissionModel()

    #Add the Isotope Contribution to the model.
    #Set the neutron speed (Ei = 0,5,14), the fraction of each
    #Isotope, and the fraction uncertainty.
    modelFullSpec.AddContribution(isotope=U235, Ei=0, fraction=1/4)
    modelFullSpec.AddContribution(isotope=U238, Ei=0, fraction=1/4)
    modelFullSpec.AddContribution(isotope=Pu240, Ei=0, fraction=1/4)
    modelFullSpec.AddContribution(isotope=Pu241, Ei=0, fraction=1/4)
    #Initialize your summation engine, and add the
    #FissionModel to your summation engine
    resultFullSpec = SumEngine()
    resultFullSpec.AddModel(modelFullSpec)

    betaSpectraDBFullSpec = BetaEngine(resultFullSpec.FPYlist.keys())
    betaSpectraDBFullSpec.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1, lower=-1.0, thresh=0.0, erange=20.0)

    resultFullSpec.CalcReactorSpectrum(betaSpectraDBFullSpec)
    print(resultFullSpec.reactorSpectrum)
    #Draw the model and save it as a png file.
    result.Draw("equal Model")


    print(resultFullSpec.missingCont)
    print(resultFullSpec.missingBranch)

#-------------------------------------------------------------------

    #Neutrino Spectrum Built from individual fission fragments-

    modelIndBB = FissionModel()

    modelIndBB.AddIstp(37,92,0.074)
    modelIndBB.AddIstp(37,93,0.037)
    modelIndBB.AddIstp(38,95,0.026)
    modelIndBB.AddIstp(39,96,0.136)
    modelIndBB.AddIstp(39,97,0.038)
    modelIndBB.AddIstp(41,100,0.03)
    modelIndBB.AddIstp(55,140,0.027)
    modelIndBB.AddIstp(55, 142, 0.05)
    resultIndBB = SumEngine()
    resultIndBB.AddModel(modelIndBB)
    betaSpectraDBIndBB = BetaEngine(resultIndBB.FPYlist.keys())
    betaSpectraDBIndBB.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1, lower=1.0, thresh=0.0, erange=20.0)

    resultIndBB.CalcReactorSpectrum(betaSpectraDBIndBB)
    print(resultIndBB.reactorSpectrum)
    resultIndBB.Draw("Dwyer-Langford Isotopes")



#-----------------------------------------------------------------

    #Neutrino Spectrum from a combination of Isotopes and individual Branches

    modelComb = FissionModel()
    modelComb.AddContribution(isotope=U235, Ei=0, Fraction= 0.5)
    modelComb.AddContribution(isotope=Pu240, Ei=0, Fraction=0.25)
    modelComb.AddIstp(39,96, 0.125)
    modelComb.AddIstp(37,92,0.0125)

    resultComb = SumEngine()
    resultComb.AddModel(modelComb)
    betaSpectraDBComb = BetaEngine(resultComb.FPYlist.keys())
    betaSpectraDBComb.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1, lower=1.0, thresh=0.0, erange=20.0)

    resultComb.CalcReactorSpectrum(betaSpectraDBComb)
    print(resultComb.reactorSpectrum)
    resultComb.Draw("Combined Model")

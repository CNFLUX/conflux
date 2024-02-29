import matplotlib.pyplot as plt
import numpy as np
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import json
import time





def pullJson(file):
    with open(file) as f:
        data = json.load(f)
    
    CalcType = data["Calculation_Type"]
    E_range = data["Energy_Range"]
    Isos = data["Isotopes"]
    fracs = data["Fraction"]
    return CalcType, E_range, Isos, fracs


def LoadIso(isotopes):
    IsoDict = {}
    print(isotopes)
    for i in isotopes:
        if(str(i) == "U235"):
            U235 = FissionIstp(92, 235)
            U235.LoadFissionDB()
            IsoDict[i] = U235
        elif(str(i) == "U238"):
            U238 = FissionIstp(92, 238)
            U238.LoadFissionDB()
            IsoDict[i] = U238
        elif(str(i) == "Pu239"):
            Pu239 = FissionIstp(94, 239)
            Pu239.LoadFissionDB()
            IsoDict[i] = Pu239
        elif(str(i) == "Pu241"):
            Pu241 = FissionIstp(94,241)
            Pu241.LoadFissionDB()
            IsoDict[i] = Pu241
        else:
            print("One of your isotopes is not" + 
                " a valid isotope. Please edit and try again")
    return IsoDict  

def LoadFisModel(isotopes, frac):
    Isodic = LoadIso(isotopes)
    FisModel = FissionModel()
    counter = 0
    for i in Isodic:
        if(i == "U235" or i == "Pu239" or i == "Pu241"):
            FisModel.AddContribution(Isodic[i], Ei=0, fraction=frac[counter])
            counter+= 1
            print("We added " + str(i) + " To our model")
        elif(i =="U238"):
            FisModel.AddContribution(Isodic[i], Ei=0.5, fraction=frac[counter])
            counter +=1
        else:
            print("You inputted a wrong isotope, please edit" + 
            " your json file and try again")
    return FisModel

def CalcSumModel(FisMod, Erange):
    result = SumEngine(xbins = np.arange(Erange[0], Erange[1], 0.1))
    result.AddModel(FisMod)
    betaSpectraDB = BetaEngine(result.FPYlist.keys(), xbins = np.arange(Erange[0], Erange[1], 0.1))
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=Erange)
    result.CalcReactorSpectrum(betaSpectraDB)
    return result

def dataOut(erange, spectra, CalcType):
    print("This is a data output script")
    now = time.time()
    text = str(CalcType) + "_" +  str(now)
    with open(text + ".txt", "w") as file:
        file.write("[E], [Spectrum]" + "\n")
        for i in range(len(spectra)):
            file.write(str(erange[i]) + "," + str(spectra[i]) + "\n")
    file.close()
    plt.plot(erange, spectra)
    plt.savefig(text + ".png")


def setCalcType(CalcType):
    if(str(CalcType) == "Summation"):
        print("Summation Mode")
    if(str(CalcType) == "Conversion"):
        print("Conversion Mode")
    if(str(CalcType) == "Direct"):
        print("Direct Calculation Mode")
    else:
        print("please check the mode you inputted.")


if __name__ == "__main__":

    Calctype, Erange, iso, frac = pullJson("testJson.json")
    
    FisModel = LoadFisModel(iso, frac)
    spectra = CalcSumModel(FisModel, Erange)
    dataOut(spectra.xbins, spectra.spectrum, Calctype)
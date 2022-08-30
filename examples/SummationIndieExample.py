import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import operator

# conflux modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine

if __name__ == "__main__":
    spectBE = [0.05, 10.05]
    width = 0.1
    model = FissionModel()
    
    inputname = sys.argv[1]
    with open(inputname) as inputfile:
        reader = csv.DictReader(inputfile, dialect='excel', delimiter=',')
        for row in reader:
            Z = int(row["Z"])
            A = int(row["A"])
            I = int(row["I"])
            Y = float(row["Y"])
            Yerr = float(row["Yerr"])
            model.AddIstp(Z, A, fraction=Y, isomer=I, d_frac = Yerr)
            
    sum1 = SumEngine(binwidths=width, spectRange=spectBE)
    sum1.AddModel(model)

    betaSpectraDB = BetaEngine(sum1.FPYlist.keys(), binwidths=width, spectRange=spectBE)
    #betaSpectraDB = BetaEngine(newlist)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])

    sum1.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False)
    summed_spect = sum1.spectrum
    summed_err = sum1.uncertainty
    summed_model_err = sum1.modelUnc
    summed_yerr = sum1.yieldUnc

    print(sum1.totalYield)
    print(sum1.missingCount)
    print(sum1.missingBranch)
    
    sum1.SaveToFile('235U_nu_jeff_indies.csv')

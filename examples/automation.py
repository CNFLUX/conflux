import sys
import numpy as np
import matplotlib.pyplot as plt

from BetaEngine import BetaEngine
from FPYEngine import FissionModel, FissionIstp
from SumEngine import SumEngine


def autobranches(fileName):
    Raw_Data = []
    f = open(fileName, 'r')
    for i in f:
        temp = i.split(";")
        Raw_Data.append(temp)
    f.close()
    for i in range(len(Raw_Data)):
        Raw_Data[i][0] = Raw_Data[i][0][6:]

    NewData = []
    runNum = 0
    for i in Raw_Data:
        IndRuns = []
        IndRuns.append(Raw_Data[runNum][0])
        runNum = runNum +1

        for j in i:
            if len(j) == 1:
                pass
            else:
                processData = j.split(",")
                length = len(processData)
                processData[0] = processData[0][1:]
                processData[length - 2] = processData[length - 2][1:]
                processData[length - 1] = processData[length - 1][:3]
                IndRuns.append(processData)
        NewData.append(IndRuns)

    for i in NewData:
        Run_Number = i[0]
        model = FissionModel()
        for j in range(1, len(i)):
            model.AddIstp(int(i[j][0]), int(i[j][1])
                          , float(i[j][2]), float(i[j][3]))
        result = SumEngine()
        result.AddModel(model)
        betaSpectra = BetaEngine(result.FPYlist.keys())
        betaSpectra.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1, thresh=0.0, erange=20.0)

        result.CalcReactorSpectrum(betaSpectra)
        result.Draw(Run_Number)
def autoIso(fileName):
    Raw_Data = []

    f = open(file, 'r')
    for i in f:
        temp = i.split(";")
        Raw_Data.append(temp)
    f.close()

    for i in range(len(Raw_Data)):
        Raw_Data[i][0] = Raw_Data[i][0][6:]
        Raw_Data[i][1] = Raw_Data[i][1][6:]

    NewData = []
    runNum = 0

    for i in Raw_Data:
        indRuns = []
        indRuns.append(Raw_Data[runNum][0])
        indRuns.append(Raw_Data[runNum][1])
        runNum = runNum + 1

        for j in i:
            if len(j) == 1:
                pass
            else:
                processData = j.split(",")
                processData[0] = processData[0][1:]
                processData[1] = processData[1][1:]
                processData[2] = processData[2][:3]
                indRuns.append(processData)
        NewData.append(indRuns)

    for i in NewData:
        Run_Number = i[0]
        model = FissionModel()
        for j in range(2, len(i)):

            if ((i[j][0] == "235" or i[j][0] == "238") and (float(i[j][1]) > 0)):
                Istp = FissionIstp(92, int(i[j][0]))
                Istp.LoadDB()
                model.AddContribution(isotope=Istp, Ei=0, fraction=float(i[j][1]), d_frac=float(i[j][2]))
            elif ((i[j][0] == "239" or i[j][0] == "241") and (float(i[j][i]) > 0)):
                Istp = FissionIstp(94, int(i[j][0]))
                Istp.LoadDB()
                model.AddContribution(isotope=Istp, Ei=0, fraction=float(i[j][1]), d_frac=float(i[j][2]))

        result = SumEngine()
        result.AddModel(model)

        betaSpectra = BetaEngine(result.FPYlist.keys())
        betaSpectra.CalcBetaSpectra(nu_spectrum=True, binwidths=0.1, lower=1.0, thresh=0.0, erange=20.0)

        result.CalcReactorSpectrum(betaSpectra)
        print(result.reactorSpectrum)
        result.Draw(str(i[0]) + ":Power- " + str(i[1]))
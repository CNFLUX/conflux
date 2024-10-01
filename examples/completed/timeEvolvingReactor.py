from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import numpy as np
import matplotlib.pyplot as plt
import csv

if __name__ == "__main__":
    #Load in the Fission Data for our fissile isotopes.
    #Also load up the energy range array
    U235 = FissionIstp(92, 235, Ei=0)
    U235.LoadFissionDB()
    U238 = FissionIstp(92, 238, Ei=0.5)
    U238.LoadFissionDB()
    Pu239 = FissionIstp(94, 239, Ei=0)
    Pu239.LoadFissionDB()
    Pu241 = FissionIstp(94, 241, Ei=0)
    Pu241.LoadFissionDB()

    e = np.arange(0,20., 0.5)

    #Load up the BetaSpectrum, and Calculate the beta spectrum for our fission products.
    betaSpectraDB = BetaEngine(xbins=e)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)
    U235.CalcBetaSpectra(betaSpectraDB)
    U238.CalcBetaSpectra(betaSpectraDB)
    Pu239.CalcBetaSpectra(betaSpectraDB)
    Pu241.CalcBetaSpectra(betaSpectraDB)


    #Open up the reactor output, set Days the reactor is running
    file = open(CONFLUX_DB+'/example_models/timeEvolvingReactor.csv')
    days = [0.0, 0.1, 2.0, 20.0, 100.0, 200.0, 300.0, 400.0, 600.0,
            800.0, 1000.0, 1200.0, 1400.0, 1600.0]

    #read out the first line which should be the
    #isotope names that are being fed into our simulation
    csvreader = csv.reader(file)
    print(next(csvreader))

    step = 0
    fig = plt.plot()
    for row in csvreader:
        #Figure out what day we are doing the calculation for
        name = days[step]
        #Initialize an array for the storage of our raw fraction data.
        raw = []
        total = 0
        #This is to skip over certain days inside our reactor, so that
        #We do simulation calculations for only specific days
        if step != 0 and step != 13 and step != 7 and step != 4 and step != 10:
            print("step is: " + str(step))
            step += 1
            continue
        percentage = []
        #This is to figure out the isotope fractions from the raw inputed data.
        for num in row:
            total += float(num)
            raw.append(float(num))
        for num in raw:
            per = num / total
            #Check for if our percentage is 0, change it to something extremely small
            if (per == 0):
                per = .000000001
            percentage.append(per)

        #Initialize the Summation Engine, Add all the fission isotopes to the engine with their
        #Fractional contributions to the total spectrum.
        SummationEngine = SumEngine(betaSpectraDB)
        SummationEngine.AddFissionIstp(U235, "U235", count = percentage[0])
        SummationEngine.AddFissionIstp(U238, "U238", count = percentage[1])
        SummationEngine.AddFissionIstp(Pu239, "Pu239", count = percentage[2])
        SummationEngine.AddFissionIstp(Pu241, "Pu241", count = percentage[3])
        SummationEngine.CalcReactorSpectrum()

        #Plot the spectrum at the current timestep
        sumX = SummationEngine.xbins
        sumY = SummationEngine.spectrum
        Labelmaker = "Reactor on for " + str(name) + " days"
        plt.plot(sumX, sumY, label = Labelmaker)
        step += 1

    #Plotting
    plt.yscale('log')
    plt.legend()
    plt.xlabel("E (in MeV)")
    plt.ylabel("neutrinos/MeV/Fission")

    plt.savefig("TimeEvolution")
    print("This is working as intended")

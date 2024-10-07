from conflux.BetaEngine import BetaEngine, CONFLUX_DB
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import numpy as np
import matplotlib.pyplot as plt
import csv

if __name__ == "__main__":
    # Load in the Fission Data for our fissile isotopes.
    # Also load up the energy range array
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

    #Create a summation engine, load in the fissile isotopes with an initial contribution of 0
    SummationEngine = SumEngine(betaSpectraDB)
    SummationEngine.AddFissionIstp(U235, "U235", count = 0)
    SummationEngine.AddFissionIstp(U238, "U238", count = 0)
    SummationEngine.AddFissionIstp(Pu239, "Pu239", count = 0)
    SummationEngine.AddFissionIstp(Pu241, "Pu241", count = 0)

    #Next, Open up the reactor output file that contains the % - fissile contribution
    file = open(CONFLUX_DB+'/example_models/timeEvolvingReactor.csv')

    #read out the first line which should be the isotope names that are being fed into our simulation
    #As well as the # of Days (This is the Header)
    csvreader = csv.reader(file)
    print(next(csvreader))

    step = 0
    fig = plt.plot()
    #Next, I iterate through the lines in the csv file
    for row in csvreader:
        #Figure out what day we are doing the calculation for
        name = float(row[0])
        #Initialize an array for the storage of our raw fraction data.
        raw = []
        total = 0
        #This is to skip over certain days inside our reactor, so that
        #We do simulation calculations for only specific days
        if step != 0 and step != 13 and step != 7 and step != 4 and step != 10:
            print("We are skipping over the information from " + str(step) + " Days after the reactor was turned on")
            step += 1
            continue
        percentage = []
        #This is to figure out the isotope fractions from the raw inputed data.
        for num in range(len(row)):
            #Skip over the first column of data, which containes the number of days after reactor on
            if (num == 0):
                continue
            else:
                total += float(row[num])
                raw.append(float(row[num]))
        for num in raw:
            per = num / total
            #Check for if our percentage is 0, change it to something extremely small
            if (per == 0):
                per = .00000000000000000000000001
            percentage.append(per)

        #Edit the fractional contribution from each fissile isotope, and then calculate the 
        #Total Reactor spectrum.
        SummationEngine.EditContribution(istpname="U235", count = percentage[0], d_count = 0)
        SummationEngine.EditContribution(istpname="U238", count = percentage[1], d_count = 0)
        SummationEngine.EditContribution(istpname="Pu239", count = percentage[2], d_count = 0)
        SummationEngine.EditContribution(istpname="Pu241", count = percentage[3], d_count = 0)
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

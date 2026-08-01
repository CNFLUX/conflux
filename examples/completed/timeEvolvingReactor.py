from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
import matplotlib.pyplot as plt
import csv

if __name__ == "__main__":
    #Load in the Fission Data
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB()
    U238 = FissionIstp(92, 238)
    U238.LoadFissionDB()
    Pu239 = FissionIstp(94, 239)
    Pu239.LoadFissionDB()
    Pu241 = FissionIstp(94, 241)
    Pu241.LoadFissionDB()

    #Open up the reactor output, set Days the reactor is running
    file = open('./timeEvolvingReactor.csv')
    days = [0.0, 0.1, 2.0, 20.0, 100.0, 200.0, 300.0, 400.0, 600.0,
            800.0, 1000.0, 1200.0, 1400.0, 1600.0]

    #read out the first line which should be the
    #isotope names that are being fed into our simulation
    csvreader = csv.reader(file)
    print(next(csvreader))

    step = 0
    fig = plt.plot()
    for row in csvreader:
        name = days[step]
        raw = []
        total = 0
        #This is to skip over certain days inside our reactor, so that
        #We do simulation calculations for only specific days
        if step != 0 and step != 13 and step != 7 and step != 4 and step != 10:
            print("step is: " + str(step))
            step += 1
            continue
        percentage = []
        #This is to figure out the isotope fractions
        #From the raw inputed data.
        for num in row:
            total += float(num)
            raw.append(float(num))
        for num in raw:
            per = num / total
            percentage.append(per)

        #This is the normal Summation calculation
        #that's been done in the Summation Engine
        FisModel = FissionModel()

        result = SumEngine()

        FisModel.AddContribution(isotope=U235, Ei = 0, fraction = percentage[0])
        FisModel.AddContribution(isotope=U238, Ei = 0.5, fraction = percentage[1])
        FisModel.AddContribution(isotope=Pu239, Ei = 0, fraction = percentage[2])
        FisModel.AddContribution(isotope=Pu241, Ei = 0, fraction = percentage[3])
        result.AddModel(FisModel)

        betaSpectraDB = BetaEngine(result.FPYlist.keys())
        betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.,20.0])
        result.CalcReactorSpectrum(betaSpectraDB)

        sumX = result.xbins
        sumY = result.spectrum
        lab = "Reactor on for " + str(name) + " days"
        plt.plot(sumX, sumY, label = lab)
        step += 1

    #Plotting
    #plt.yscale('log')
    plt.legend()
    plt.xlabel("E (in MeV)")
    plt.ylabel("neutrinos/MeV/Fission")

    plt.savefig("TimeEvolution")
    print("This is working as intended")
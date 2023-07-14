from conflux.FPYEngine import FissionIstp, FissionModel
from conflux.BetaEngine import BetaEngine

import csv

# test on how to define fission isotopes and add them to the reactor model
if __name__ == "__main__":
    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB()
    U235.LoadCovariance()
    # U235.LoadCorrelation('/Users/zhang39//Downloads/')
    # U235.CalcCovariance(Ei =0)
    with open('cov_235_U_processed_v2.csv', 'w', newline='') as csvoutput:
        fieldnames = ['']
        for i in U235.CFPY[0]:
            fieldnames.append(i)
        writer = csv.DictWriter(csvoutput, fieldnames=fieldnames)
        writer.writeheader()
        for i in U235.CFPY[0]:
            rowcontent = U235.CFPY[0][i].cov
            rowcontent[''] = i
            writer.writerow(rowcontent)


    Pu239 = FissionIstp(94, 239)
    Pu239.LoadFissionDB()

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0, fraction=1, d_frac=0.0)
    #model.AddContribution(isotope=Pu239, Ei = 0, fraction=0.4, d_frac=0.05)
    #model.AddIstp(390960)
    totalyield = 0
    betaSpectraDB = BetaEngine(model.FPYlist.keys())
    betaSpectraDB.LoadBetaDB()
    for ZAI in betaSpectraDB.istplist:
        print("ZAI", ZAI)
    totalFPY = 0
    for FPZAI in model.FPYlist:
        if FPZAI in betaSpectraDB.istplist and betaSpectraDB.istplist[FPZAI].missing == False:
            print('nuclide: ', FPZAI, 'y: ', model.FPYlist[FPZAI].y, 'yerr: ', model.FPYlist[FPZAI].yerr )
            totalyield += model.FPYlist[FPZAI].y
            totalFPY += 1
    print("totalyield", totalyield)
    print("totalFPY", totalFPY)

    model.DrawBranches("frac.4.png")
    model.SaveToFile('frac.4.csv')

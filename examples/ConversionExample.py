# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData

# test
if __name__ == "__main__":
    beta235 = BetaData("./data/conversionDB/U_235_e_2014.csv")
    beta239 = BetaData("./data/conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData("./data/conversionDB/Pu_241_e_2014.csv")
    #print(beta235.y)
    #print(beta235.yerr)

    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)
    U235.LoadFissionDB()
    Pu239.LoadFissionDB()
    Pu241.LoadFissionDB()

    # vbtest = VirtualBranch(U235)
    # vbtest.CalcZAavg(6,7)
    # print(vbtest.Aavg, vbtest.Zavg)

    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)
    convertmodel.VBfit(0.5)
    convertmodel.DrawVB("U235", "U235_conversion.png")

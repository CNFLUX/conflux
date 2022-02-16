# local modules
from conflux.BetaEngine import BetaEngine, BetaBranch
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.ConversionEngine import ConversionEngine

# test
if __name__ == "__main__":
    beta235 = BetaData("./conversionDB/U_235_e_2014.csv")
    beta239 = BetaData("./conversionDB/Pu_239_e_2014.csv")
    beta241 = BetaData("./conversionDB/Pu_241_e_2014.csv")
    #print(beta235.y)
    #print(beta235.yerr)

    U235 = FissionIstp(92, 235)
    Pu239 = FissionIstp(94, 239)
    Pu241 = FissionIstp(94, 241)
    U235.LoadDB()
    Pu239.LoadDB()
    Pu241.LoadDB()

    # vbtest = VirtualBranch(U235)
    # vbtest.CalcZAavg(6,7)
    # print(vbtest.Aavg, vbtest.Zavg)

    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, Pu239, "Pu239", 1.0)
    convertmodel.VBfit(0.5)
    convertmodel.DrawVB("Pu239", "Pu239_convert_test_0.5.png")

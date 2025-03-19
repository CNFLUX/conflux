from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp, FissionModel



if __name__ == "__main__":
    # Fission Model is a mostly defunct class. It is not required when running the Summation
    # Or conversion calculations. However, it is useful when storing reactor informatoin
    # Of both Fission and non-fission Isotopes. Let's see how it works.

    #To start, I'm going to initialize two Fission isotopes, U235 and Pu239, and load in their
    #Correlation Fission information, with an incident neutron energy of 0
    U235 = FissionIstp(92,235,Ei=0)
    U235.LoadFissionDB()
    U235.LoadCorrelation()
    Pu239 = FissionIstp(94, 239, Ei=0)
    Pu239.LoadFissionDB()
    Pu239.LoadCorrelation()

    # Next, I'm going to go ahead and put both of these classes into a Fission Model. First,
    # I initialize my fission model, and then I add in both Fission Isotopes. Note,
    # I need to specify the contribution of each of these Fission Isotopes to the overall
    # model.
    FisModel = FissionModel()
    FisModel.AddContribution(U235, fraction=0.46, d_frac=0.001)
    FisModel.AddContribution(Pu239, fraction=0.46, d_frac=0.001)

    # I can also go ahead and add some other random non-beta decaying isotopes into the model.
    # Let's add C-12 as well, but a very small amount without any uncertainty in the amount that's
    # in there. 
    FisModel.AddIstp(Z = 6, A = 12, fraction=.08)

    #Finally, let me go ahead and draw the beta branches of the resultant model, as well as save.
    #All the FPY information to a file. Note, I didn't actually have to calculate the spectrum of
    #Anything in order to do this.
    FisModel.DrawBranches("Branches.png")
    FisModel.SaveToFile("FPYs.csv")
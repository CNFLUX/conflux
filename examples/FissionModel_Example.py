from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp, FissionModel



if __name__ == "__main__":
    #So, I am going to assume that you've looked at FissionIstp_Example
    #Before coming to this example, as the first chunk of code that I'm going
    #to run is going to be identical to the code in that example, and so if you
    #Would like a deep dive as to why I'm running specific commands, or why I'm
    #calling the commands in the order that I am, I suggest you take a look there.

    #------Identical to FissionIstp_Example-----#
    U235 = FissionIstp(92,235)
    U235.LoadFissionDB()
    U235.LoadCovariance()
    #-------------------------------------------#

    #Okay, now that that's out of the way, let's get onto the new part of this.
    #First, we're going to initialize a Fissionmodel, and then we're going to add
    #The contribution from U235 to our fission model. Next, we're going to
    #Draw the individual beta branches of U235. Finally, we're going to get
    #Some information about an isotope based off its ZAI number (Refer to the
    #BetaEngine_Example to understand what the ZAI number is.


    #Here I'm initializing my model
    model = FissionModel()

    #And here I'm adding the U235 Isotope that I had initialized before into the
    #Fission Model. Couple of things I should mention here. The first parameter of this
    #is the isotope, the second is the neutron energy, third is the fractional contribution
    #of that specific isotope on the overall reactor, and finally the fourth parameter
    #Is the fractional uncertainty of your isotope. Note, you can omit the last variable if
    #the uncertainty is unknown. The code will still behave itself.
    #model.AddContribution(isotope=U235, Ei=0, fraction=1.0, d_frac=0.0)
    model.AddContribution(U235, 0, fraction=1.0, d_frac=0.0)

    #Here, I'm drawing all the individual beta branches for U235. As of right now, there's
    #no customizability to the draw functions, though that might change in future updates.
    model.DrawBranches("U235_Beta_Branch.png")

    #Finally, we are going to get the nuclide information for Bromine-90. Like I said earlier,
    #If there's any confusion as to what the ZAI number of an element is, please refer back to
    #The BetaEngine_Example function.

    isotope = model.GetNuclide(350900)
    print("This is my atomic number:", isotope.Z," And my mass number: ", isotope.A)
    print("This is the ZAI number: ", isotope.FPZAI)
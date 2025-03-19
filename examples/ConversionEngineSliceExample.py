from conflux.FPYEngine import CONFLUX_DB, FissionIstp
from conflux.ConversionEngine import ConversionEngine, BetaData
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    #Here, I want to show off CONFLUXs' ability to define the beta branch size when
    #Carrying out the conversion calculation. The size of the beta branch determines
    #How well we fit our beta spectrum before converting to a neutrino spectrum.
    #Most of the code here is the same as the ConversionEngine_Example.py file,
    #So I will not be going in depth into the parts that are the same. For more questions,
    #Please look at that file.

    #Load up the energy range for the calculation (0 to 10 MeV with 100keV steps)
    e = np.linspace(0., 10., 100)

    #Load the beta information from the ConversionDB directory
    #(change this to the local directory for U_235_e_2014.csv)
    beta235 = BetaData(CONFLUX_DB+"/conversionDB/U_235_e_2014.csv")

    #Initialize a U235 Fission Isotope with incident neutron energy of 0 (Thermal)
    U235 = FissionIstp(92,235, Ei=0)
    U235.LoadFissionDB()

    #Load up the conversion model, add the beta data as well as the fission isotpe to
    #The model. the fractional contribution of U235 in this case is 1.
    convertmodel = ConversionEngine()
    convertmodel.AddBetaData(beta235, U235, "U235", 1.0)

    #Run the virtual beta branch fitting for incrementing branch slices.
    #Branch slices are increased by 0.2, up to 4.0
    fig = plt.figure()
    for slicesize in [0.25, 0.5, 1]:
        print(slicesize)
        #This is some code for setting up our plots.        
        # plt.yscale('log')
        plt.xlabel("E (MeV)")
        plt.ylabel("neutrino/fission/MeV")
        #This is the size of our virtual branch slices
        convertmodel.VBfitbeta("U235", slicesize=slicesize)
        #The rest of the plotting is the same as what you see in the ConversionEngine_Example.py file.
        name = f"U235_slice_size_{slicesize}_MeV"
        plt.errorbar(e, convertmodel.vblist["U235"].SumBranches(e, nu_spectrum = True), label=f"slice size {slicesize}_MeV")
        plt.legend()
        
        # OPTIONAL below uncertainty calculation

        #This calculates the errors of our conversion calcualtion, with 50 samples.
        #I generate a covariance matrix for both the betas and the neutrinos
        # covmat = convertmodel.vblist["U235"].Covariance(beta235, e, nu_spectrum=False, samples=50)
        # covmat_nu = convertmodel.vblist["U235"].Covariance(beta235, e, nu_spectrum=True, samples=50)

        #Here, I calculate the beta spectrum, and then calculate the relative errors of it from the
        #Covariance matrices I generated above.
        # spectrumBeta = convertmodel.vblist["U235"].SumBranches(e, nu_spectrum=False)
        # relativeErr = np.zeros(len(e))
        # for i in range(len(spectrumBeta)):
        #     if spectrumBeta[i] > 0:
        #         relativeErr[i] = np.sqrt(covmat[i][i]) / spectrumBeta[i]

        #Here, I calculate the neutrino spectrum, and then calculate the relative errors of it from the
        #Covariance matrices I generated above.
        # spectrumNu = convertmodel.vblist["U235"].SumBranches(e, nu_spectrum=True)
        # relativeErrNu = np.zeros(len(e))
        # for i in range(len(spectrumNu)):
        #     if spectrumNu[i] > 0:
        #         relativeErrNu[i] = np.sqrt(covmat_nu[i][i]) / spectrumNu[i]


        #Finally, the rest of this is plotting the relative errors for both the betas and the neutrinos for our given
        #Energy range
        # plt.ylim([-.3, .3])
        # plt.fill_between(e, relativeErr, -relativeErr, label='beta', alpha=0.4)
        # plt.fill_between(e, relativeErrNu, -relativeErrNu, label='neutrino', alpha=0.4)
        # plt.legend()
        
    name = "U235_slice.pdf"
    fig.savefig(name)

    plt.close()
import numpy as np
import matplotlib.pyplot as plt

# local modules
from conflux.BetaPlusEngine import BetaEngine

from xraydb import XrayDB
from scipy.constants import alpha

xdb = XrayDB()

eV_2_MeV = 1e-6

# def EC_nu_spectrum(xbins, istp):
#     Z = istp.Z
#     Zd = Z - 1 # the atom number of Z 
#     xray_edges = xdb.xray_edges(Z)
#     # neutrino energy to be filled
#     E_nu = []
#     decay_frac = np.array([])
#     decay_frac_unc = np.array([])
#     # branch fraction of all neutrino energy
    
#     shellnum = {"K":1, "L":2, "M":3, "N":4, "O":5, "P":6}
#     suborbit_hi = [2, 3, 4, 5]
#     for E0, branch in istp.branches.items():
#         ecfrac = branch.ecfrac
#         sigma_ecfrac = branch.sigma_ecfrac
#         print(ecfrac, sigma_ecfrac)
#         shellfrac = []
#         for shell, xray_edge in xray_edges.items():
#             omega = 0
#             Eb = xray_edge.energy*eV_2_MeV
#             Knu = E0-Eb 
            
#             # exclude forbidden capture where E0 is smaller than the binding energy.
#             if Knu < 0: 
#                 continue 

#             # identify the subshell where electron is captured
#             orbitnum = 1
#             hi_suborbit = 1
#             for key, item in shellnum.items():
#                 if key in shell: orbitnum = item
#             for i in suborbit_hi:
#                 if str(i) in shell: hi_suborbit = i
            
#             # calculate the effective fraction ofthe 
#             if hi_suborbit == 1:
#                 omega = Zd**3/orbitnum**3
#             elif hi_suborbit == 2:
#                 omega = Zd**5/2**5 * alpha**2            

#             # add the fraciton and neutrino energy of the subshell capture
#             shellfrac.append(omega)
#             E_nu.append(E0 - Eb)
        
#         # normalize the fraction
#         shellfrac = np.array(shellfrac)
#         decay_frac = np.append(decay_frac, shellfrac/sum(shellfrac)*ecfrac)
#         decay_frac_unc = np.append(decay_frac_unc, shellfrac/sum(shellfrac)*sigma_ecfrac)
    
#     binwidth = xbins[1]-xbins[0]
#     y = np.zeros(len(xbins))
#     y_unc = np.zeros(len(xbins))
#     for i in range(len(E_nu)):
#         E0 = E_nu[i]
#         I0 = decay_frac[i]
#         sigmaI = decay_frac_unc[i]
#         if not (0.0 <= E0 <= xbins[-1]): 
#             continue
#         t = (E0 - xbins[0]) / binwidth                     # position in bins
#         i = int(np.floor(t))
#         frac = t - i
#         if i < 0:
#             y[0] += I0
#             y_unc[0] += sigmaI
#         elif i >= y.size - 1:
#             y[-1] += I0
#             y_unc[-1] += sigmaI
#         else:
#             y[i]     += I0 * (1.0 - frac)
#             y[i + 1] += I0 * frac
#             y_unc[i]     += sigmaI * (1.0 - frac)
#             y_unc[i + 1] += sigmaI * frac
    
#     y = y/binwidth
#     y_unc = y_unc/binwidth
#     return y, y_unc

if __name__ == "__main__":
    #One thing I didn't specify in the BetaEngine_Example is that you can in fact input a 
    #Custom list of isotopes whose Beta Spectra you want to calculate. You do not need to 
    #Load and Calculate every single Fission Product, and for certain calculations, it's
    #Less resource intensive to input a list.

    #The reaction starts in much the same way as it did in the BetaEngine_Example.py file, 
    #Initialize some energy range, load up a BetaEngine, and then calculate the BetaSpectra
    #However, before I Initialize the BetaEngine, I will define a list of isotopes by their
    #ZAI number. These isotopes will be the ones whose betaSpectrum I want to calculate.

    #initialize energy
    e = np.arange(0.,3.,0.01)

    #Load up a list of isotopes I want to calculate
    istplist = [
        40070, 
        60100,
        90180, 
        220450,
        521180, 
        822000
        ]

    #Load up the BetaEngine, and pass the isotope list above into the BetaEngine
    BetaSpectraDB = BetaEngine(istplist, xbins = e)

    #Next, I will go ahead and Calculate the BetaSpectrum, and plot out the result, 
    #As well as the total spectrum
    BetaSpectraDB.CalcBetaSpectra(nu_spectrum=True)

    #Iterate over the isotopic list inside the BetaEngine we intitialized, and plot out
    #the Beta Spectrum for each isotope.
    fig = plt.plot()
    for i in istplist:
        betaIstp = BetaSpectraDB.istplist[i]
        plt.errorbar(e, betaIstp.spectrum, betaIstp.uncertainty, label = betaIstp.name)

    plt.xlabel("E (MeV)")
    plt.ylabel("neutrino/MeV")
    plt.legend()
    # plt.show()
    plt.savefig("Individual_bp_neutrino_spectrum.pdf")
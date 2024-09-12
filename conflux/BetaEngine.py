import numpy as np
from scipy import special, interpolate, integrate
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm

from conflux.config import CONFLUX_DB
from conflux.Basic import Spectrum, Summed
from conflux.bsg.Constants import ELECTRON_MASS_MEV, NATURAL_LENGTH
from conflux.bsg.SpectralFunctions import (phase_space, 
                                           fermi_function, 
                                           finite_size_L0, 
                                           recoil_gamow_teller, 
                                           radiative_correction, 
                                           radiative_correction_neutrino,
                                           recoil_Coulomb_gamow_teller,                               
                                           atomic_screening,
                                           shape_factor_gamow_teller,
                                           shape_factor_unique_forbidden)
from conflux.bsg.Functions import getEltonNuclearRadius
from conflux.bsg.Screening import screening_potential

#########################################
# Final neutrino and antineutrino spectra

def electron(ebeta, p, numass=0):
    result = 0.
    W = ebeta/ELECTRON_MASS_MEV + 1
    result = (phase_space(W, numass=numass, **p, )
            *fermi_function(W, **p)
            *finite_size_L0(W, **p)
            *recoil_gamow_teller(W, **p)
            *radiative_correction(W, **p)
            *recoil_Coulomb_gamow_teller(W, **p)
            *atomic_screening(W, **p)
            )

    if p['L'] == 0:
        result *= shape_factor_gamow_teller(W, **p)
    else:
        result *= shape_factor_unique_forbidden(W, **p)

    return result

def neutrino(enu, p, numass=0):
    result = 0.
    W0 = p['W0']
    Wv = W0-enu/(ELECTRON_MASS_MEV*1.0) #enu/ELECTRON_MASS_MEV + 1
    result = (phase_space(Wv, numass=numass, **p)
            *fermi_function(Wv, **p)
            *finite_size_L0(Wv, **p)
            *recoil_gamow_teller(Wv, **p)
            *radiative_correction_neutrino(Wv, **p)
            *recoil_Coulomb_gamow_teller(Wv, **p)
            *atomic_screening(Wv, **p)
            )

    if p['L'] == 0:
        result *= shape_factor_gamow_teller(Wv, **p)
    else:
        result *= shape_factor_unique_forbidden(Wv, **p)
    return result

# BetaBranch class to save the isotopic information
class BetaBranch(Spectrum):
    def __init__(self, Z, A, I, Q, E0, sigma_E0, frac, sigma_frac,
                forbiddenness=0, bAc=4.7, xbins=np.arange(0, 20, 0.1),
                custom_func=None):
        self.ID = E0

        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q
        self.ZAI=Z*1e4+A*10+I

        self.xbins = xbins
        self.nbin = len(xbins)
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

        self.E0 = E0
        self.sigma_E0 = sigma_E0
        self.frac = frac
        self.sigma_frac = sigma_frac

        self.forbiddenness = forbiddenness
        
        self.custom_func = custom_func
        self.Parameters = {
            'Z': Z,
            'A': A,
            'W0': E0/ELECTRON_MASS_MEV + 1,
            'R': getEltonNuclearRadius(A) * 1e-15 / NATURAL_LENGTH,
            'L': forbiddenness,
            'c': 1.0,
            'b': bAc*A,
            'd': 0.0,
            'Lambda': 0.0,
            'l': screening_potential(Z)
        }

        self.corr = {E0:1}  # correlation with other branches of the same isotope
        self.cov = {E0:self.sigma_frac**2}

    # display the vital info of branch
    def Display(self):
        print("Branch E0 = " + str(self.E0)+ "+\-" + str(self.sigma_E0)
            + ", frac = " + str(self.frac) + "+\-"+str(self.sigma_frac))

    # set correlation of this branch fraction and all another branches
    def SetCovariance(self, otherBranch, correlation):
        if self.E0 == otherBranch.E0:
            return
        self.corr[otherBranch.E0] = correlation
        self.cov[otherBranch.E0] = (self.sigma_frac * correlation
                                    * otherBranch.sigma_frac)

    # beta spectrum shape as function of energy
    def BetaSpectrum(self, x, nu_spectrum=False, numass=0):
        Parameters = deepcopy(self.Parameters)

        # prevent out-of-range (> Q value)variable to create insane results
        rangeCorrect = x <= self.E0

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, Parameters, numass=numass)
        else:
            function = lambda x: electron(x, Parameters, numass=numass)
            
        if (self.custom_func!=None):
            function = lambda x: self.custom_func(x, Parameters, numass=numass)

        result = function(x)
        result = np.nan_to_num(result, nan=0.0)
        return result*rangeCorrect

    # calculate the spectrum uncertainty with MC sampling
    def SpectUncertMC(self, x, nu_spectrum=False, samples = 30):
        if self.sigma_E0 == 0:
            return 0
        E0range = np.random.normal(self.E0, self.sigma_E0, samples)
        newParameters = deepcopy(self.Parameters)
        newParameters['W0'] = E0range/ELECTRON_MASS_MEV+1

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, newParameters)
        else:
            function = lambda x: electron(x, newParameters)
            
        if (self.custom_func!=None):
            function = lambda x: self.custom_func(x, newParameters)

        fE0 = function(x)
        fE0 = np.nan_to_num(fE0, nan=0.0)
        # Commenting out unknown function
        # if (fE0.all() < 1e-8):
        #     return 0
        return np.std(fE0)

    # bined beta spectrum
    def BinnedSpectrum(self, nu_spectrum=False):
        # to prevent lower limit of the energy range < 0
        lower = self.xbins[0]
        if (lower > self.E0):
            return 1
        # TODO make the code compatible to uneven binning
        binwidths = self.xbins[1]-self.xbins[0]
        # integrating each bin
        for k in range(0, self.nbin):
            x_low = lower
            x_high = lower+binwidths
            if x_high > self.E0:
                x_high = self.E0
            thisbinwidth = abs(x_high-x_low)
            relativewidth = abs(x_high-x_low)/binwidths
            self.spectrum[k] = (thisbinwidth*relativewidth
                                *self.BetaSpectrum(x_low, nu_spectrum))
            self.uncertainty[k] = (thisbinwidth*relativewidth
                                *self.SpectUncertMC(x_low, nu_spectrum))
            if x_high == self.E0:
                break
            lower += binwidths

        full_range = np.arange(0, 20, 0.01)
        this_range = np.arange(self.xbins[0], self.xbins[-1], 0.01)
        full_spect = self.BetaSpectrum(full_range, nu_spectrum)
        this_spect = self.BetaSpectrum(this_range, nu_spectrum)

        # normalizing the spectrum
        norm = (self.spectrum.sum()*full_spect.sum()/this_spect.sum() 
                if self.E0 > binwidths else self.spectrum.sum())

        if self.spectrum.sum() <=0:
            self.spectrum = np.zeros(self.nbin)
        else:
            self.spectrum /= norm*binwidths
            self.uncertainty /= norm*binwidths

        return 0

# class to save isotope information, including Z A I Q and beta brances
class BetaIstp(Spectrum, Summed):
    def __init__(self, Z, A, I, Q, HL, name, xbins=np.arange(0, 20, 0.1)):
        self.ZAI=Z*1e4+A*10+I # unique ID of a isotope
        self.ID = self.ZAI
        self.HL = HL

        self.xbins = xbins
        self.nbin = len(xbins)
        self.spectrum = np.zeros(self.nbin)
        self.uncertainty = np.zeros(self.nbin)

        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q
        self.name = name
        self.missing=False
        self.branches = {}

    def AddBranch(self, branch):
        """
        Add beta branch to this isotope
        Parameters:
            branch(BetaBranch)
        Returns:
            None
        """
        self.branches[branch.ID] = branch


    def EditBranch(self, defaultE0, **kwargs):
        """
        Add or edit branches to the isotope with analyzer's assumptions
        Parameters:
            E0 (float): 
                existing end point energy of the missing branch
            newE0 (float): 
                redefined end point energy of the missing branch
            fraction (float): 
                assumed fraction of the edited branch
            sigma_E0 (float): 
                1-sigma error of the input E0, default as 0
            sigma_frac (float): 
                1-sigma error of the input franction, default as 0
        Returns:
            None
        """
        # if sigma_E0 > E0:
        #     sigma_E0 = E0

        # setting up default values
        sigma_E0 = 0
        fraction = 1
        sigma_frac = 1
        forbiddenness = 0
        bAc = 4.7

        for key, value in kwargs.items():
            if key == 'E0':
                defaultE0 = value
            if key == 'sigma_E0':
                sigma_E0 = value
            if key == 'fraction':
                fraction = value
            if key == 'forbiddenness':
                forbiddenness = value
            if key == 'bAc':
                bAc = value
            if key == 'custom_func':
                bAc = custom_func

            if defaultE0 not in self.branches.keys():
                self.branches[defaultE0] = BetaBranch(self.Z, 
                                                      self.A, 
                                                      self.I, 
                                                      self.Q, 
                                                      defaultE0, 
                                                      sigma_E0, 
                                                      fraction, 
                                                      sigma_frac,
                                                      forbiddenness,
                                                      bAc=bAc, 
                                                      xbins=self.xbins,
                                                      custom_func=custom_func)

            if hasattr(self.branches[defaultE0], key):
                setattr(self.branches[defaultE0], key, value)

    def MaxBranch(self):
        """
        Get the maximum E0 branch of the isotope
        Returns:
            self.branch(Emax) (BetaBranch)
        """
        Emax = max(self.branches)
        return self.branches[Emax]

    def CalcCovariance(self, GSF=True):
        """
        Calculate the covaraince matrix for all the beta branches of this
        isotope
        Parameters:
            GSF (boolean): determine whether to calculate covaraince matrix 
            with ground state feeding
        Returns:
            None
        """

        MaxBranch = self.MaxBranch()
        # obtain the ground state branch info
        totalFrac = 0
        GSFrac = MaxBranch.frac
        for E0, branch in self.branches.items():
            totalFrac += branch.frac
        restFrac = totalFrac-GSFrac

        if GSF == False or MaxBranch.E0 != self.Q or restFrac == 0:
            for i, branchi in self.branches.items():
                for j, branchj in self.branches.items():
                    branchi.SetCovariance(branchj, 0)
            return

        # calculate the covariance matrix with ground state anticorrelation
        for i, branchi in self.branches.items():
            for j, branchj in self.branches.items():
                if i == self.Q:
                    correlation = -1*branchj.frac/restFrac
                    branchi.SetCovariance(branchj, correlation)
                    branchj.SetCovariance(branchi, correlation)
                elif j != self.Q:
                    branchi.SetCovariance(branchj, 0)

    def SumSpectra(self, nu_spectrum=True, branchErange=[0, 20.]):
        """
        Calculate the cumulative beta/antineutrino spectrum of all branches
        Parameters:
            nu_spectrum (boolean): 
                Determines if you are calculating a beta spectra or a neutrino 
                spectra
            binwidths (float): 
                The width of the bins used in creating your spectra
            spectRange (list): 
                Define the upper and power bounds of the spectrum
        Returns:
            None
        """

        self.spectrum=np.zeros(self.nbin)
        self.uncertainty=np.zeros(self.nbin)
        self.spectUnc=np.zeros(self.nbin) # theoretical uncertainty
        self.branchUnc=np.zeros(self.nbin)
        self.totalUnc=np.zeros(self.nbin)

        # calculate the total uncertaintty and append spectra
        for E0i, branchi in self.branches.items():
            if(E0i < branchErange[0] or E0i > branchErange[1]):
                continue
            si = branchi.spectrum
            fi = branchi.frac
            # di = branchi.sigma_frac

            branchi.BinnedSpectrum(nu_spectrum)
            self.spectrum += si*fi
            self.spectUnc += branchi.uncertainty*fi

            for E0j, branchj in self.branches.items():
                sj = branchj.spectrum
                # fj = branchi.frac
                # dj = branchj.sigma_frac
                cov_bij = branchi.cov[E0j]
                sigma_bij = si*cov_bij*sj
                self.branchUnc += sigma_bij
                if (E0i==E0j):
                    self.totalUnc += (branchi.uncertainty*fi)**2 + sigma_bij
                else:
                    self.totalUnc += sigma_bij


        self.branchUnc = np.sqrt(self.branchUnc)
        self.totalUnc = np.sqrt(self.totalUnc)
        self.uncertainty = self.totalUnc

    def Display(self):
        """
        Display isotope property and branch information
        """
        print('Beta isotope: '+self.name+', ZAI = '+str(self.ZAI)+', Q = '
            +str(self.Q)+" MeV, " +str(len(self.branches))+" branches")
        for E0, branch in self.branches.items():
            branch.Display()

    def decay_time_adjust(self, begin=0, end=0):
        '''
        function to calculate the decay rate of this isotope from the begin of
        counting to the end of the counting.
        if no argument is given, return 1
        '''
        # print(self.name, self.HL)
        if begin < end:
            return 2**(-begin/self.HL)-2**(-end/self.HL)
        else:
            return 1

# BetaEngine tallys beta branches in the betaDB and calculate theoretical beta 
# spectra of all tallied branches
# if inputlist is not given, load the entire betaDB from the default betaDB
class BetaEngine:
    def __init__(self, 
                 inputlist=None, 
                 targetDB=CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml",
                 xbins=np.arange(0, 20, 0.1),
                 custom_func=None):
        self.inputlist = inputlist
        self.istplist = {}
        self.xbins = xbins

        self.LoadBetaDB(targetDB)   # loadBetaDB automatically
        
        self.custom_func=custom_func

    def LoadBetaDB(self, targetDB=CONFLUX_DB+"/betaDB/ENSDFbetaDB2.xml"):
        """Load default or input betaDB to obtain beta decay informtion
        """
        useInputList = True # test if the engine is defined with an inputlist
        if self.inputlist == None:
            print("Loading all beta data from the default betaDB...")
            self.inputlist = []
            useInputList = False

        print(targetDB)
        print("Searching DB: "+targetDB+"...")
        print("Loading spectra of beta branches...")

        tree = ET.parse(targetDB)
        root = tree.getroot()
        for isotope in root:
            ZAI = int(isotope.attrib['isotope'])
            Q = float(isotope.attrib['Q'])
            HL = float(isotope.attrib['HL'])
            name = isotope.attrib['name']

            # if input list is not given, include all isotopes
            if not useInputList:
                self.inputlist.append(ZAI)

            if (ZAI in self.inputlist):
                Z = int(ZAI/10000)
                A = int(ZAI%10000/10)
                I = int(ZAI%10)

                # if same ZAI reappear it is tagged as an isomer (unmarked in
                # the original DB)
                if ZAI in self.istplist:
                    I += 1
                    ZAI += 1

                betaIstp = BetaIstp(Z, A, I, Q, HL, name, xbins=self.xbins)

                # Adding missing branches below
                if len(isotope) < 1:
                    betaIstp.missing = True
                    betaIstp.EditBranch(betaIstp.Q, E0=betaIstp.Q, fraction=1)
                    self.istplist[ZAI] = betaIstp
                    continue

                # some isotopes contain summed branch fraction greater than 1
                # prepare to normalize
                fracsum = 0
                for branch in isotope:
                    fraction = float(branch.attrib['fraction'])
                    fracsum += fraction

                # actually assign values from database to branches
                for branch in isotope:
                    E0 = float(branch.attrib['end_point_E'])
                    if E0<=0:
                        continue
                    sigma_E0 = float(branch.attrib['sigma_E0'])
                    spin_par_changes = (branch.attrib['dJpi'])

                    # converting spin change to forbidden types
                    ftypes = [['0', '1'], ['0-', '1-', '2-'], ['2', '3'],
                              ['3-', '4-'], ['4', '5']]
                    firstftypes = [-10, -11, 10]
                    forbiddenness = 6
                    for i in range(len(ftypes)):
                        for j in range(len(spin_par_changes)):
                            if (spin_par_changes[j] in ftypes[i] and 
                                i < abs(forbiddenness)):
                                forbiddenness = -i
                                if i == 1:
                                    forbiddenness = firstftypes[j]
                                elif spin_par_changes[j] == ftypes[i][-1]:
                                    forbiddenness = i

                    # assign fraction values to branches
                    # normalize if greater than one
                    fraction = float(branch.attrib['fraction'])
                    sigma_frac = float(branch.attrib['sigma_frac'])
                    if fracsum > 1:
                        fraction /= fracsum
                        sigma_frac /= fracsum

                    betaBranch = BetaBranch(Z, A, I, Q, E0, sigma_E0, fraction,
                                            sigma_frac, forbiddenness,
                                            xbins=self.xbins
                                            custom_func=self.custom_func)
                    betaIstp.AddBranch(betaBranch)

                if betaIstp.branches:
                    self.istplist[ZAI] = betaIstp


    def CalcBetaSpectra(self, targetDB = None, nu_spectrum=True,
                        branchErange=[0.0, 20.0], GSF=False, silent=False):
        """
        Calculates beta spectra of the list of beta-decaying isotopes

        Parameters:
            targetDB (String): 
                the path to the betaSpectra database. if none, use the default database.
            nu_spectrum (boolean): 
                Determines if you are calculating a beta spectra or a neutrino spectra
            binwidths (float): 
                The width of the bins used in creating your spectra.
            spectRange (list): 
                define the upper and power bounds of the spectrum
            branchErange (list): 
                defind the range of interested Q values
        Returns:
            None

        """

        for ZAI in tqdm(self.istplist, desc="Calculating beta/neutrino spectra of "+str(len(self.istplist))+ " isotopes", disable=silent):
            betaIstp = self.istplist[ZAI]
            if betaIstp.Q < branchErange[0] or betaIstp.Q > branchErange[1]:
                continue

            betaIstp.CalcCovariance(GSF=GSF)

            betaIstp.SumSpectra(nu_spectrum, branchErange)

if __name__ == "__main__":
    x = np.arange(0, 10, 0.1)
    binwidth = 1

    testlist = [390960, 390961, 521331, 531371, 922390, 932390]
    testEngine = BetaEngine(xbins=x)
    testEngine.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20], GSF=False)

    y2 = testEngine.istplist[390960].spectrum
    y2err = testEngine.istplist[390960].uncertainty
    plt.figure()
    plt.xlabel("E (keV)")
    plt.errorbar(x/1e-3, y2, label=str(390960)+"_nu", yerr=y2err)
    plt.legend()
    plt.show()

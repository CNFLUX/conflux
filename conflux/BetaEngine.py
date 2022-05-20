import sys
import os
import argparse
import numpy as np
from scipy import constants, special, interpolate, integrate
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from copy import deepcopy
import pkg_resources
import timeit

######################
# Declaring constants

ELECTRON_MASS_MEV = constants.m_e/constants.eV*constants.speed_of_light**2/1E6
ELECTRON_MASS_EV = constants.m_e/constants.eV*constants.speed_of_light**2
PROTON_MASS_EV = constants.m_p/constants.eV*constants.speed_of_light**2
NUCLEON_MASS_EV = (constants.m_p+constants.m_n)/2.0/constants.eV*constants.speed_of_light**2
NUCLEON_MASS_W = NUCLEON_MASS_EV/(1.0*ELECTRON_MASS_EV)
PROTON_MASS_W = PROTON_MASS_EV/(1.0*ELECTRON_MASS_EV)
FERMI_to_W= 1e-15/(constants.hbar*constants.speed_of_light/constants.e)*ELECTRON_MASS_EV

################
# Fermi theory
def nuclear_radius(A):
    result = 1.121*pow(A,1/3.0)+2.426*pow(A,-1/3.0)-6.614/A
    return result

WO = lambda energy: energy/(1.0*ELECTRON_MASS_MEV) + 1.0

p = lambda W: np.sqrt(W**2 - 1.)

y = lambda W, Z: 1.0*Z*W/p(W)*constants.fine_structure

gamma = lambda Z: np.sqrt(1.0 - (constants.fine_structure*Z*1.0)**2)

phasespace = lambda W, W0: p(W)*W*(W0-W)*(W0-W)

def F(y, gamma, p, R):
    result = 0.
    res = special.loggamma(gamma+y*1j)
    absgamma2= np.exp(2*res.real)
    g2g = special.gamma(2*gamma+1)
    result = 2*(gamma+1)*absgamma2*np.exp(y*np.pi)/(pow(2*p*R,(2*(1-gamma)))*g2g**2)
    return result

###############
# Finite size

l0dat = [[0.115, -1.8123, 8.2498, -11.223, -14.854, 32.086],
		    [-0.00062, 0.007165, 0.01841, -0.53736, 1.2691, -1.5467],
		    [0.02482, -0.5975, 4.84199, -15.3374, 23.9774, -12.6534],
		    [-0.14038, 3.64953,-38.8143,172.1368,-346.708, 288.7873],
		    [0.008152,-1.15664,49.9663,-273.711,657.6292, -603.7033],
		    [1.2145, -23.9931,149.9718,-471.2985, 662.1909,-305.6804],
		    [-1.5632,33.4192,-255.1333,938.5297,-1641.2845,1095.358]]

def L0(W, Z, R, gamma):
    result = 0.
    alpha = constants.fine_structure
    result = 1. + (13/60.)*(alpha*Z)*(alpha*Z) - (W*R*alpha*Z*(41 - 26*gamma))/(15*(2*gamma - 1)) - (alpha*Z*R*gamma*(17 - 2*gamma))/(30*W*(2*gamma -1)) - 0.41*(R-0.0164)*pow(alpha*Z,4.5);
    return result

def af(j0, Z):
    j = j0+1
    result = 0.
    for i in range(0, 6):
        result += l0dat[j][i]*pow(Z*constants.fine_structure, i+1.0)
    return result

def L0b(W, Z, R):
    result = 0.
    for k in range(0, 6):
        result += af(k,Z)*pow(W*R,k)
    result=result+af(-1,Z)*R/W;
    return result

######################
# Screening Correction
#TODO: build a screening effect database to replace

Wb = lambda energy, V0: WO(energy)-V0

nn_energy =[1,8,13,16,23,27,29,49,84,92]
nn_mfp=[1,1.42,1.484,1.497,1.52,1.544,1.561,1.637,1.838,1.907]
nn_spline = interpolate.InterpolatedUnivariateSpline(nn_energy,nn_mfp)

def NN(Z):
    return nn_spline(Z)

def V0(Z):
    return NN(Z-1)*constants.fine_structure**2*pow(Z-1, 4/3.)

def S(energy, Z):
    result = 0.
    result = (Wb(energy,V0(Z))/WO(energy))*pow(p(Wb(energy,V0(Z)))/p(WO(energy)),2*gamma(Z)-1)*np.exp(-np.pi*(y(WO(energy),Z)-y(Wb(energy,V0(Z)),Z)))
    res = special.loggamma((gamma(Z) + 1j*y(Wb(energy,V0(Z)),Z)))
    absgamma2b= np.exp(2*res.real)
    res = special.loggamma((gamma(Z) + 1j*y(WO(energy),Z)))
    absgamma2= np.exp(2*res.real)
    result = result*absgamma2b/absgamma2
    return result

#####################
# Shape Corrections

def CC(R, Z, W, W0):
    result = 0.
    alpha = constants.fine_structure*1.0
    result += (-(233/630.0))*(alpha*Z)*(alpha*Z) - (1/5.0)*(W0*R)*(W0*R) + (2.0/35)*W0*R*alpha*Z
    result += W*((-(21/35.0))*R*alpha*Z + (4*W0*R*R)/9.0)
    result += -((4*R*R)/9.0)*W*W
    return result+1

###########################
# QED radiative corrections

def G(W, W0):
    result = 0.
    beta = p(W)/W
    result += 3*np.log(NUCLEON_MASS_W)-3/4.0+4*(np.arctanh(beta)/beta-1.0)*((W0-W)/(3*W)-3/2.0+np.log(2*(W0-W)))
    result += (4.0*special.spence(1-(2*beta)/(1+ beta)))/beta
    result += (np.arctanh(beta)*(2.0*(1.0+beta**2)+(W0-W)*(W0-W)/(6.0*W**2)-4.0*np.arctanh(beta)))/beta
    result = result*constants.fine_structure/(2.0*np.pi)+1.0
    return result

def GN(W):
    result=0.
    beta=np.sqrt(W*W-1.0)/W
    result += 3.0*np.log(PROTON_MASS_W)+23.0/4.0-8.0/beta * special.spence(1-2.0 * beta/(1+beta))
    result += 8.0*(np.arctanh(beta)/beta-1.0)*np.log(2.0*W*beta)
    result += 4.0*np.arctanh(beta)/beta*((7.0+3.0*beta*beta)/8.0-2.0*np.arctanh(beta))
    result = result*constants.fine_structure/(2.0*np.pi)+1.0
    return result

#######################
# Forbidden shapes

def forbidden(W, W0, WM, ftype):
    p1 = p(W)
    pn = W0-W
    result = 1.0
    enu = pn*ELECTRON_MASS_MEV*1.0
    ee = W*ELECTRON_MASS_MEV*1.0
    beta = p1/W
    shape = 0.
    wm = 0.
    pe = np.sqrt(ee**2 - ELECTRON_MASS_MEV**2*1.0)

    if ((ftype==1) or (ftype==-2)): # first unique, 2nd non-unique
        result = (pn*pn+p1*p1)*(1+W*ELECTRON_MASS_MEV*WM)
    if ((ftype==2) or (ftype==-3)): # 2nd unique, 3rd non-unique
        result=(pow(pn,4)+10.0/3*pn*pn*p1*p1+pow(p1,4))*(1+W*ELECTRON_MASS_MEV*WM)
    if ((ftype==3) or (ftype==-4)): # 3rd unique, 4th non-unique
        result=(pow(pn,6)+7.0*pow(pn,4)*p1*p1+7.0*pow(p1,4)*pn*pn+pow(p1,6))*(1+W*ELECTRON_MASS_MEV*WM)
    if (ftype == -10):  #   first non-unique 0-
        shape=(pe*pe+enu*enu+2*beta*beta*enu*ee)
        wm=0
        wm=wm*WM
        result=(1+wm)*shape
    if (ftype==-11):    #   first non-unique 1-
        shape= pe*pe + enu*enu - 4.0/3.0*beta*beta*enu*ee
        wm=((pe*pe + enu*enu)*(beta*beta*ee - enu) + 2.0*beta*beta*ee*enu*(enu - ee)/3.0)/((pe*pe + enu*enu - 4.0*beta*beta*enu*ee/3.0))
        wm=wm*WM
        result=(1+wm)*shape
    if (ftype == 10):   #   first unique, 2-
        shape=pe**2 + enu**2
        wm=3.0/5.0*((pe*pe + enu*enu)*(beta*beta*ee - enu) + 2.0*beta*beta*ee*enu*(enu - ee)/3.0)/((pe*pe + enu*enu))
        wm=wm*WM
        result=(1+wm)*shape

    return result

#########################################
# Final neutrino and antineutrino spectra

class Parameters_t:
    def __init__(self, e0, Z, A, ftype, WM, thresh=0):
        self.e0 = e0
        self.Z = Z
        self.A = A
        self.ftype = ftype
        self.WM = WM
        self.thresh = thresh

def neutrino(enu, Parameters):
    np.seterr(divide='ignore')
    buf = Parameters
    result = 0.
    R= 1.0*FERMI_to_W * nuclear_radius(buf.A)
    W0=WO(buf.e0)
    W=W0-enu/(ELECTRON_MASS_MEV*1.0)
    thresh=(W0-1.0-V0(buf.Z))*(ELECTRON_MASS_MEV*1.0)
    buf.thresh = thresh
    result = forbidden(W,W0,buf.WM,buf.ftype)*GN(W)*(L0(W, buf.Z, R, gamma(buf.Z))+L0b(W,buf.Z,R))*CC(R, buf.Z, W, W0)*phasespace(W, W0)*F(y(W,buf.Z), gamma(buf.Z), p(W), R)

    belowThresh = enu<thresh   # prepared for a list of comparison
    #print("E value: ",buf.e0, "calc: ", result, "stdev: ", np.std(result))
    return result*(S(buf.e0-enu, buf.Z)**belowThresh)

def electron(ebeta, Parameters):
    buf = Parameters
    result = 0.
    R = 1.0*FERMI_to_W * nuclear_radius(buf.A)
    W0=WO(buf.e0)
    W=WO(ebeta)
    thresh=V0(buf.Z)*(ELECTRON_MASS_MEV*1.0)
    buf.thresh = thresh
    result = forbidden(W,W0,buf.WM,buf.ftype)*G(W,W0)*(L0(W, buf.Z, R, gamma(buf.Z))+L0b(W,buf.Z,R))* CC(R, buf.Z, W, W0)*phasespace(W, W0)*F(y(W,buf.Z), gamma(buf.Z), p(W), R)
    
    aboveThresh = (ebeta>=thresh)   # prepared for a list of comparison
    return result*(S(ebeta,buf.Z)**aboveThresh)

def integral(nu_spectrum, p, x_low, x_high):
    erg = 0.0
    xh = 0.0
    x_low = max([x_low, 1e-6])

    if (x_low >= p.e0-1e-6):
        return 0
    if (x_high > p.e0):
        xh = p.e0
    else:
        xh = x_high

    if (nu_spectrum == True):
        function = lambda x: neutrino(x, p)
    else:
        function = lambda x: electron(x, p)

    #function(x_low)
    if (x_low <= p.thresh and xh >= p.thresh): # wierd error because of integrater when difference is smaller than a value

        ergl = [0.0, 0.0]
        ergh = [0.0, 0.0]

        # a quick integration mathod through trapezoid algorithm
        mid = (x_low+p.thresh)/2
        ergl = abs(p.thresh-x_low)*function(mid)
        mid = (xh+p.thresh)/2
        ergh = abs(p.thresh-xh)*function(mid)
        erg = ergl+ergh

        return erg

    mid = (xh+x_low)/2
    erg = abs(xh-x_low)*function(mid)
    
    return erg

    #print(x_low, p.thresh, x_high, erg[0])

# BetaBranch class to save the isotopic information
class BetaBranch:
    def __init__(self, Z, A, I, Q, E0, sigma_E0, frac, sigma_frac, forbiddeness=0, WM=0.0047):
        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q
        self.ZAI=Z*1e4+A*10+I

        self.E0 = E0
        self.sigma_E0 = sigma_E0
        self.frac = frac
        self.sigma_frac = sigma_frac

        self.forbiddeness = forbiddeness
        self.WM = WM
        
        self.Parameters = Parameters_t(Z = self.Z+1, A = self.A, e0=self.E0, WM=self.WM, ftype=self.forbiddeness)
        
        self.corr = {E0:1}  # correlation with other branches of the same isotope
        self.cov = {E0:self.sigma_frac**2}
        
    # display the vital info of branch
    def Display(self):
        print("Branch E0 = " +str(self.E0)+"+\-"+str(self.sigma_E0)+", frac = "+str(self.frac)+"+\-"+str(self.sigma_frac))
    
    # set correlation of this branch fraction and all another branches
    def SetCovariance(self, otherBranch, correlation):
        if self.E0 == otherBranch.E0:
            return
        self.corr[otherBranch.E0] = correlation
        self.cov[otherBranch.E0] = self.sigma_frac*correlation*otherBranch.sigma_frac
            
    # beta spectrum shape as function of energy
    def BetaSpectrum(self, x, nu_spectrum=False):
        Parameters = deepcopy(self.Parameters)
        rangeCorrect = x <= self.E0 # prevent out-of-range variable to output insane results

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, Parameters)
        else:
            function = lambda x: electron(x, Parameters)

        result = function(x)
        result = np.nan_to_num(result, nan=0.0)
        return result*rangeCorrect

    # calculate the spectrum uncertainty
    def SpectUncert(self, x, nu_spectrum=False):
        
        numbers = 5
        E0range = np.linspace(self.E0-self.sigma_E0, self.E0+self.sigma_E0, numbers)
        newParameters = deepcopy(self.Parameters)
        newParameters.e0 = E0range

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, newParameters)
        else:
            function = lambda x: electron(x, newParameters)
            
        fE0 = function(x)
        fE0 = np.nan_to_num(fE0, nan=0.0)
        # Commenting out unknown function
        # if (fE0.all() < 1e-8):
        #     return 0

        grad_E0 = np.gradient(fE0)
        
        return grad_E0[int(numbers/2)]*self.sigma_E0

    # calculate the spectrum uncertainty with MC sampling
    def SpectUncertMC(self, x, nu_spectrum=False, samples = 50):
        
        E0range = np.random.normal(self.E0, self.sigma_E0, samples)
        newParameters = deepcopy(self.Parameters)
        newParameters.e0 = E0range
        
        if (nu_spectrum == True):
            function = lambda x: neutrino(x, newParameters)
        else:
            function = lambda x: electron(x, newParameters)
            
        fE0 = function(x)
        fE0 = np.nan_to_num(fE0, nan=0.0)
        # Commenting out unknown function
        # if (fE0.all() < 1e-8):
        #     return 0
        return np.std(fE0)
        
    # bined beta spectrum
    def BinnedSpectrum(self, nu_spectrum=False, binwidths=0.1, spectRange=[-1.0, 20.0]):
        bins = int(spectRange[1]/binwidths)
        self.result = np.zeros(bins)
        self.uncertainty = np.zeros(bins)

        lower = spectRange[0]
        if (lower > self.E0):
            return 1
        if (lower<0):
            lower=binwidths/2.0

        # integrating each bin
        # startTiming = timeit.default_timer()

        for k in range(0, bins):
            x_low = lower
            x_high = lower+binwidths
            if x_high > self.E0:
                x_high = self.E0
            self.result[k] = abs(x_high-x_low)*self.BetaSpectrum((x_low+x_high)/2, nu_spectrum)
            #self.uncertainty[k] = abs(x_high-x_low)*(self.SpectUncert(x_low, nu_spectrum)+self.SpectUncert(x_high, nu_spectrum))/2
            #gradUnc = self.uncertainty[k]
            #print("gradUnc", gradUnc)
            self.uncertainty[k] = abs(x_high-x_low)*self.SpectUncertMC((x_low+x_high)/2, nu_spectrum)
            if x_high == self.E0:
                break
            #MCUnc = self.uncertainty[k]
            #print("MCUnc", MCUnc)
            lower+=binwidths

        # endTiming = timeit.default_timer()
        # runTime = endTiming-startTiming
        # print("runtime", runTime)
        
        # normalizing the spectrum
        norm = self.result.sum()
        #print(self.result, self.uncertainty, norm*binwidths)
        if norm <=0:
            self.result =np.zeros(bins)
        else:
            self.result /= norm*binwidths
            self.uncertainty /= norm*binwidths
            #print(self.result, self.uncertainty, norm*binwidths)
        return 0

# class to save isotope information, including Z A I Q and beta brances
class BetaIstp:
    def __init__(self, Z, A, I, Q, name):
        self.Z = Z
        self.A = A
        self.I = I
        self.Q = Q
        self.name = name
        self.ZAI=Z*1e4+A*10+I
        self.branch={}
        self.missing=False
    
    def AddBranch(self, branch):
        """
        Add beta branch to this isotope
        Returns:
            None
        """
        self.branch[branch.E0] = branch
    
    def EditBranch(self, E0, fraction, sigma_E0 = 0., sigma_frac = 0., forbiddeness = 0):
        """
        Add or edit branches to the isotope with analyzer's assumptions
        Parameters:
            E0 (float): assumed end point energy of the missing branch
            fraction (float): assumed fraction of the edited branch
            sigma_E0 (float): 1-sigma error of the input E0, default as 0
            sigma_frac (float): 1-sigma error of the input franction, default as 0
        Returns:
            None
        """
        self.branch[E0] = BetaBranch(self.Z, self.A, self.I, self.Q, E0, sigma_E0, fraction, sigma_frac, forbiddeness)
        
    def MaxBranch(self):
        """
        Get the maximum E0 branch of the isotope
        Returns:
            self.branch(Emax) (BetaBranch)
        """
        Emax = max(self.branch)
        return self.branch[Emax]
    
    def CalcCovariance(self, GSF=True):
        """
        Calculate the covaraince matrix for all the beta branches of this isotope
        Parameters:
            GSF (boolean): determine whether to calculate covaraince matrix with ground state feeding
        Returns:
            None
        """
        MaxBranch = self.MaxBranch()
        # obtain the ground state branch info
        totalFrac = 0
        GSFrac = MaxBranch.frac
        for E0, branch in self.branch.items():
            totalFrac += branch.frac
        restFrac = totalFrac-GSFrac
        
        if GSF == False or MaxBranch.E0 != self.Q or restFrac == 0:
            for i, branchi in self.branch.items():
                for j, branchj in self.branch.items():
                    branchi.SetCovariance(branchj, 0)
            return
        
        # calculate the covariance matrix with ground state anticorrelation
        for i, branchi in self.branch.items():
            for j, branchj in self.branch.items():
                if i == self.Q:
                    correlation = -1*branchj.frac/restFrac
                    branchi.SetCovariance(branchj, correlation)
                    branchj.SetCovariance(branchi, correlation)
                elif j != self.Q:
                    branchi.SetCovariance(branchj, 0)
                
    def CalcBetaSpectrum(self, nu_spectrum=True, binwidths=0.1, spectRange=[-1.0, 20.0]):
        """
        Calculate the cumulative beta/antineutrino spectrum of all branches
        Parameters:
            nu_spectrum (boolean): Determines if you are calculating a beta spectra or a neutrino spectra
            binwidths (float): The width of the bins used in creating your spectra
            spectRange (list): define the upper and power bounds of the spectrum
        Returns:
            None
        """
        spectLow = spectRange[0] if spectRange[0] > 0 else 0
        spectHigh = spectRange[1]
        bins = int((spectHigh-spectLow)/binwidths)
        self.spectrum=np.zeros(bins)
        self.spectUnc=np.zeros(bins)
        self.branchUnc=np.zeros(bins)
        self.totalUnc=np.zeros(bins)
        
        for E0,branch in self.branch.items():
            branch.BinnedSpectrum(nu_spectrum, binwidths, spectRange)
            self.spectrum += branch.result*branch.frac
            self.spectUnc += branch.uncertainty*branch.frac
        
        for E0i, branchi in self.branch.items():
            si = branchi.result
            di = branchi.sigma_frac
            for E0j, branchj in self.branch.items():
                sj = branchj.result
                dj = branchj.sigma_frac
                cov_bij = branchi.cov[E0j]
                sigmab_ij = si*cov_bij*sj
                self.branchUnc += sigmab_ij
                if (E0i==E0j):
                    self.totalUnc += (branchi.uncertainty*branchi.frac)**2 + sigmab_ij
                else:
                    self.totalUnc += sigmab_ij
        self.branchUnc = np.sqrt(self.branchUnc)
        self.totalUnc = np.sqrt(self.totalUnc)
        
    def Display(self):
        """
        Display isotope property and branch information
        """
        print('Beta isotope: '+self.name+', ZAI = '+str(self.ZAI)+', Q = '+str(self.Q)+" MeV, " +str(len(self.branch))+" branches")
        for E0, branch in self.branch.items():
            branch.Display()
    
# BetaEngine tallys beta branches in the betaDB and calculate theoretical beta spectra
# of all tallied branches
# if inputlist is not given, load the entire betaDB from the default betaDB
class BetaEngine:
    def __init__(self, inputlist=None, targetDB=None):
        self.inputlist = inputlist
        self.defaultDB = os.environ["CONFLUX_DB"]+"/betaDB/ENSDFbetaDB.xml"
        
        self.LoadBetaDB(targetDB)   # loadBetaDB automatically
        
    def LoadBetaDB(self, targetDB=None):
        useInputList = True         # test if the engine is defined with an inputlist
        if self.inputlist == None:
            print("Loading all beta data from the default betaDB...")
            self.inputlist = []
            useInputList = False
            
        if targetDB == None:
            targetDB = self.defaultDB
        print("Searching DB: "+targetDB+"...")
        print("Loading spectra of beta branches...")

        self.istplist = {}
        
        tree = ET.parse(targetDB)
        root = tree.getroot()
        for isotope in root:
            ZAI = int(isotope.attrib['isotope'])
            Q = float(isotope.attrib['Q'])
            name = isotope.attrib['name']

            if not useInputList:
                self.inputlist.append(ZAI)
                
            if (ZAI in self.inputlist):
                #print(str(ZA)+"...")
                Z = int(ZAI/10000)
                A = int(ZAI%10000/10)
                I = int(ZAI%10)
                
                betaIstp = BetaIstp(Z, A, I, Q, name)
                
                # Adding missing branches below
                if len(isotope) < 1:
                    betaIstp.missing = True
                    betaIstp.EditBranch(betaIstp.Q, 1)
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
                    sigma_E0 = float(branch.attrib['sigma_E0'])
                    spin_par_changes = (branch.attrib['dJpi'])

                    # converting spin change to forbidden types
                    ftypes = [['0', '1'], ['0-', '1-', '2-'], ['2', '3'], ['3-', '4-'], ['4', '5']]
                    firstftypes = [-10, -11, 10]
                    forbiddeness = 1e3
                    for i in range(len(ftypes)):
                        for j in range(len(spin_par_changes)):
                            if spin_par_changes[j] in ftypes[i] and i < abs(forbiddeness):
                                forbiddeness = -i
                                if i == 1:
                                    forbideness = firstftypes[j]
                                elif spin_par_changes[j] == ftypes[i][-1]:
                                    forbideness = i
                    
                    # assign fraction values to branches
                    # normalize if greater than one
                    fraction = float(branch.attrib['fraction'])
                    sigma_frac = float(branch.attrib['sigma_frac'])
                    if fracsum > 1:
                        fraction /= fracsum
                        sigma_frac /= fracsum

                    betaBranch = BetaBranch(Z, A, I, Q, E0, sigma_E0, fraction, sigma_frac, forbiddeness)
                    betaIstp.AddBranch(betaBranch)
                self.istplist[ZAI] = betaIstp

        
    def CalcBetaSpectra(self, targetDB = None, nu_spectrum=True, binwidths=0.1, spectRange=[-1.0, 20.0], branchErange=[-1.0, 20.0]):
        """
        Calculates beta spectra of the list of beta-decaying isotopes

        Parameters:
            targetDB (String): the path to the betaSpectra database. if none, use the default database.
            nu_spectrum (boolean): Determines if you are calculating a beta spectra or a neutrino spectra
            binwidths (float): The width of the bins used in creating your spectra.
            spectRange (list): define the upper and power bounds of the spectrum
            branchErange (list): defind the range of interested Q values
        Returns:
            None
                
        """
        spectLow = spectRange[0] if spectRange[0] > 0 else 0
        spectHigh = spectRange[1]
        self.bins = np.arange(spectLow, spectHigh, binwidths)

        startTiming = timeit.default_timer()
        istpCount = 0
        for ZAI in self.istplist:
            betaIstp = self.istplist[ZAI]
            if betaIstp.Q < branchErange[0] or betaIstp.Q > branchErange[1]:
                continue

            betaIstp.CalcCovariance()
            betaIstp.CalcBetaSpectrum(nu_spectrum, binwidths, spectRange)
            istpCount += 1

        endTiming = timeit.default_timer()
        nBranch = istpCount
        runTime = endTiming-startTiming
        print("Finished calculating beta spectra of "+ str(nBranch) + " isotopes.")
        print("Processing time: "+str(runTime)+" seconds")

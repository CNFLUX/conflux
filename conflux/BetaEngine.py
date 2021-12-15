import sys
import argparse
import numpy as np
from scipy import constants, special, interpolate, integrate
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from copy import deepcopy

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

    if ((ftype==1) or (ftype==-2)):
        result = (pn*pn+p1*p1)*(1+W*ELECTRON_MASS_MEV*WM)
    if ((ftype==2) or (ftype==-3)):
        result=(pow(pn,4)+10.0/3*pn*pn*p1*p1+pow(p1,4))*(1+W*ELECTRON_MASS_MEV*WM)
    if ((ftype==3) or (ftype==-4)):
        result=(pow(pn,6)+7.0*pow(pn,4)*p1*p1+7.0*pow(p1,4)*pn*pn+pow(p1,6))*(1+W*ELECTRON_MASS_MEV*WM)
    if (ftype == -10):
        shape=(pe*pe+enu*enu+2*beta*beta*enu*ee)
        wm=0
        wm=wm*WM
        result=(1+wm)*shape
    if (ftype==-11):
        shape= pe*pe + enu*enu - 4.0/3.0*beta*beta*enu*ee
        wm=((pe*pe + enu*enu)*(beta*beta*ee - enu) + 2.0*beta*beta*ee*enu*(enu - ee)/3.0)/((pe*pe + enu*enu - 4.0*beta*beta*enu*ee/3.0))
        wm=wm*WM
        result=(1+wm)*shape
    if (ftype == 10):
        shape=pe**2 + enu**2
        wm=3.0/5.0*((pe*pe + enu*enu)*(beta*beta*ee - enu) + 2.0*beta*beta*ee*enu*(enu - ee)/3.0)/((pe*pe + enu*enu))
        wm=wm*WM
        result=(1+wm)*shape

    return result

#########################################
# Final neutrino and antineutrino spectra

class params_t:
    def __init__(self, e0, Z, A, ftype, WM, thresh=0):
        self.e0 = e0
        self.Z = Z
        self.A = A
        self.ftype = ftype
        self.WM = WM
        self.thresh = thresh

def neutrino(enu, params):
    buf = params
    result = 0.
    R= 1.0*FERMI_to_W * nuclear_radius(buf.A)
    W0=WO(buf.e0)
    W=W0-enu/(ELECTRON_MASS_MEV*1.0)
    thresh=(W0-1.0-V0(buf.Z))*(ELECTRON_MASS_MEV*1.0)
    buf.thresh = thresh
    result = forbidden(W,W0,buf.WM,buf.ftype)*GN(W)*(L0(W, buf.Z, R, gamma(buf.Z))+L0b(W,buf.Z,R))*CC(R, buf.Z, W, W0)*phasespace(W, W0)*F(y(W,buf.Z), gamma(buf.Z), p(W), R)

    if (enu<thresh):
        return result*S(buf.e0-enu, buf.Z)
    return result

def electron(ebeta, params):
    buf = params
    result = 0.
    R = 1.0*FERMI_to_W * nuclear_radius(buf.A)
    W0=WO(buf.e0)
    W=WO(ebeta)
    thresh=V0(buf.Z)*(ELECTRON_MASS_MEV*1.0)
    buf.thresh = thresh
    result = forbidden(W,W0,buf.WM,buf.ftype)*G(W,W0)*(L0(W, buf.Z, R, gamma(buf.Z))+L0b(W,buf.Z,R))* CC(R, buf.Z, W, W0)*phasespace(W, W0)*F(y(W,buf.Z), gamma(buf.Z), p(W), R)
    if(ebeta>=thresh):
        return result*S(ebeta,buf.Z)
    return result

def integral(nu_spectrum, p, x_low, x_high, verb):
    erg = 0.0
    error = 0.0
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
        test = p.thresh

        # a quick integration mathod through trapezoid algorithm
        if (verb == "quick"):
            mid = (x_low+p.thresh)/2
            ergl = abs(p.thresh-x_low)*function(mid)
            mid = (xh+p.thresh)/2
            ergh = abs(p.thresh-xh)*function(mid)
            erg = ergl+ergh
        else:
            if (x_low <= p.thresh-5e-4 ): ergl = integrate.quad(function, x_low, p.thresh) # set lower limit to 1e-6 to avoid 'nan' error
            if (xh >= p.thresh+5e-4): ergh = integrate.quad(function, p.thresh, xh)
            erg=ergl[0]+ergh[0]

        return erg

    if (verb == "quick"):
        mid = (xh+x_low)/2
        erg = abs(xh-x_low)*function(mid)
        return erg
    else:
        erg = integrate.quad(function, x_low, xh)
        return erg[0]
    #print(x_low, p.thresh, x_high, erg[0])


# BetaBranch class to save the isotopic information
class BetaBranch:
    def __init__(self, Z, A, frac, I, E0, sigma_E0, forbiddeness=0, WM=0.0047):
        self.Z = Z
        self.A = A
        self.I = I
        self.E0 = E0
        self.sigma_E0 = 0.05*E0 #sigma_E0
        self.frac = frac

        self.forbiddeness = forbiddeness
        self.WM = WM

        self.params = params_t(Z = self.Z+1, A = self.A, e0=self.E0, WM=self.WM, ftype=self.forbiddeness)

    # beta spectrum shape as function of energy
    def BetaSpectrum(self, x, nu_spectrum=False):
        params = deepcopy(self.params)
        rangecorrection = x <= self.E0 # prevent out-of-range variable to output insane results

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, params)
        else:
            function = lambda x: electron(x, params)

        result = function(x)
        result = np.nan_to_num(result, nan=0.0)
        return result

    # calculate the spectrum uncertainty
    def SpectUncert(self, x, nu_spectrum = False):
        numbers = 5
        E0range = np.linspace(self.E0-self.sigma_E0, self.E0+self.sigma_E0, numbers)
        newparams = deepcopy(self.params)
        newparams.e0 = E0range

        if (nu_spectrum == True):
            function = lambda x: neutrino(x, newparams)
        else:
            function = lambda x: electron(x, newparams)

        fE0 = function(x)
        fE0 = np.nan_to_num(fE0, nan=0.0)
        if (fE0.all() < 1e-8):
            return 0

        grad_E0 = np.gradient(fE0)

        return grad_E0[int(numbers/2)]*self.sigma_E0

    # bined beta spectrum
    def BinnedSpectrum(self, nu_spectrum=False, binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0):
        bins = int(erange/binwidths)
        self.result = np.zeros(bins)
        self.uncertainty = np.zeros(bins)

        if (lower > self.E0):
            return 1
        if (lower<0):
            lower=binwidths/2.0

        self.params.thresh = thresh

        # integrating each bin
        for k in range(0, bins):
            x_low = lower
            x_high = lower+binwidths
            if x_high > self.E0:
                x_high = self.E0
            self.result[k] = integral(True, self.params, x_low, x_high, "quick") #abs(x_high-x_low)*(self.BetaSpectrum(x_low)+self.BetaSpectrum(x_high))/2
            self.uncertainty[k] = abs(x_high-x_low)*(self.SpectUncert(x_low)+self.SpectUncert(x_high))/2
            lower+=binwidths

        norm = self.result.sum()
        if norm <=0:
            self.result =np.zeros(bins)
        else:
            self.result /= norm*binwidths
            self.uncertainty /= norm*binwidths
        return 0

# BetaEngine tallys beta branches in the betaDB and calculate theoretical beta spectra
# of all tallied branches
class BetaEngine:
    def __init__(self, inputlist, DBname ='betaDB/betaDB.xml'):
        self.inputlist = inputlist
        self.istplist = {}
        self.spectralist = {}
        self.uncertaintylist = {}
        self.DBname = DBname

    def LoadBetaDB(self):
        print("Searching DB: "+self.DBname+"...")
        print("Loading spectra of beta branches:")

        tree = ET.parse(self.DBname)
        root = tree.getroot()
        for isotope in root:
            ZAI = int(isotope.attrib['isotope'])
            if (ZAI in self.inputlist):
                #print(str(ZA)+"...")
                Z = int(ZAI/10000)
                A = int(ZAI%10000/10)
                I = int(ZAI%10)
                betaIstp = {}
                for branch in isotope:
                    E0 = float(branch.attrib['end_point_E'])
                    sigma_E0 = float(branch.attrib['sigma_E0'])
                    forbiddeness = int(branch.attrib['forbideness'])
                    fraction = float(branch.attrib['fraction'])
                    sigma_frac = float(branch.attrib['sigma_frac'])

                    betaIstp[E0] = BetaBranch(Z, A, fraction, I, E0, sigma_E0, forbiddeness)
                self.istplist[ZAI] = betaIstp

    def CalcBetaSpectra(self, nu_spectrum=True, binwidths=0.1, lower=-1.0, thresh=0.0, erange = 20.0):
        self.LoadBetaDB()
        bins = int(erange/binwidths)
        for ZAI in self.istplist:
            branchspectrum = np.zeros(bins)
            branchuncertainty = np.zeros(bins)
            for E0, branch in self.istplist[ZAI].items():
                branch.BinnedSpectrum(nu_spectrum, binwidths, lower, thresh, erange)
                branch.result *= branch.frac
                branch.uncertainty *= branch.frac

                branchspectrum += branch.result
                branchuncertainty += branch.uncertainty

            self.spectralist[ZAI] = branchspectrum
            self.uncertaintylist[ZAI] = branchuncertainty


if __name__ == "__main__":
    testlist = [521340, 531340, 350900]
    testEngine = BetaEngine(testlist)
    testEngine.CalcBetaSpectra(nu_spectrum=True)
    #print(testEngine.spectralist[521340])
    print(testEngine.spectralist[350900])
    print(testEngine.uncertaintylist[350900])

    fig = plt.figure()
    x = np.linspace(0, 20, 200)
    plt.errorbar(x, testEngine.spectralist[350900], yerr=testEngine.uncertaintylist[350900])
    #plt.draw()
    fig.savefig("errorbartest.png")

    #testbeta = BetaBranch(1, 3, 1.0, 0, )

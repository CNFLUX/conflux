#! /bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from math import *
import bisect

# conflux modules
from conflux.BetaEngine import BetaEngine, init_libBSG
from conflux.FPYEngine import FissionIstp
from conflux.IBDXSection import ibd_xsection_cm2

# fission energies from arXiv:hep-ph/0410100 Table 1
J_per_MeV = 1.60218e-13
E_fission_235 = 202.79 * J_per_MeV          # energy per fission, joules
E_fission_U9 = 207.32 * J_per_MeV           # energy per fission, joules
E_TNT = 4.184e9                             # energy per ton TNT, joules
fissions_per_ton_235 = E_TNT/E_fission_235  # number of fissions per ton TNT
fissions_per_ton_U9 = E_TNT/E_fission_U9

IBD_E0 = 1.29333236 + 0.511


# 1000 ton TNT equivalent fission burst
m_TNT = 1000

# optional use of faster compiled libBSG:
# init_libBSG()

#############
# Calculation energy grid

#xbins = np.arange(0, 13, 0.1)

#lnE0 = log(0.01)
#lnE1 = log(13.)
#xbins = [0,] + [exp(u) for u in np.arange(lnE0, lnE1, (lnE1 - lnE0)/50.)]

k = 0.05
dx0 = 0.02
xbins = [(0.2 - dx0)*2/pi*(x*atan(k*x) - log(1+(k*x)**2)/(2*k)) + dx0*x for x in range(80)]

# include IBD start energy as calculation grid point
bisect.insort_left(xbins, IBD_E0)
xbins = np.array(xbins)

betaSpectraDB = BetaEngine()
betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, uncertainties = False)

############################
############################

# Load the independent fission product yields from the JEFF database
U235 = FissionIstp(92, 235, 0.4, DB='JEFF', IFPY=True)
U235.LoadFissionDB(DB = 'JEFF')
U235.LoadCorrelation(DB = 'JEFF')

# Load the independent fission product yields from the JEFF database
P9 = FissionIstp(94, 239, 0.4, DB='JEFF', IFPY=True)
P9.LoadFissionDB(DB = 'JEFF')
P9.LoadCorrelation(DB = 'JEFF')

class TimeSpectra:
    def __init__(self):
        self.spect_time = []
        self.spect_sum_t = []
        self.spect_sum_t_ibd = []

        self.cumuspect = []
        self.cumuIBD = []

def time_spectra(fissIsot, windows, fissions = 1):

    TS = TimeSpectra()
    TS.windows = windows

    yields = fissIsot.GetYields()
    decays = betaSpectraDB.CalcDecays(yields)

    TS.totalspect = betaSpectraDB.SumBranches(decays, xbins) * fissions # the full spectrum summed with all fission products
    TS.totalspect_ibd = TS.totalspect * ibd_xsection_cm2(xbins) # the IBD detector spectrum

    # generate neutrino spectra in different time window after ignition
    for i in range(len(windows) - 1):

        # calculate the beta spectra within different windows after the fission
        decays = betaSpectraDB.CalcDecays(yields, tbegin=windows[i], tend=windows[i+1])
        spect = betaSpectraDB.SumBranches(decays, xbins) * fissions

        # add spectrum of each window to a cumulative spectrum
        if TS.cumuspect: TS.cumuspect.append(TS.cumuspect[-1] + spect)
        else: TS.cumuspect = [spect,]

        TS.cumuIBD.append(TS.cumuspect[-1] * ibd_xsection_cm2(xbins))

        TS.spect_sum_t.append(100 * sum(TS.cumuspect[-1]) / sum(TS.totalspect))
        TS.spect_sum_t_ibd.append(100 * sum(TS.cumuIBD) / sum(TS.totalspect_ibd))

    return TS

windows = [0, 1, 10, 100, 1000, 10000]
TSU = time_spectra(U235, windows,  m_TNT * fissions_per_ton_235)
TSP = time_spectra(P9, windows,  m_TNT * fissions_per_ton_U9)

fig, ax = plt.subplots()
ax.set(xlabel='E (MeV)', ylabel='neutrinos/MeV from %g ton yield'%m_TNT)
for n,C in enumerate(TSU.cumuspect):
    l0, = ax.plot(xbins, C, label="%g s"%TSU.windows[n+1])
    ax.plot(xbins, TSP.cumuspect[n], linestyle="--", color=l0.get_color())

l0, = ax.plot(xbins, TSU.totalspect, label="total")
ax.plot(xbins, TSP.totalspect, linestyle="--", color=l0.get_color())

ax.legend()
plt.ylim(0., None)
plt.xlim(0., 8)
fig.savefig("Neutrinos_time.pdf")


# IBD target parameters
rho = 1.09 # g/cm^3
nH = 74.8 * 6.02e20 # H / cm^3
H_per_kg = 1e3 * nH/rho
cxnorm = H_per_kg / (4*pi) / 100**2
dx0 = IBD_E0 - 2*0.511

fig, ax = plt.subplots()
ytitle = 'IBD/MeV / (kg scintillator / m$^2$) from %g ton yield'%m_TNT
ax.set(xlabel='IBD E (MeV), positron + annihilation', ylabel=ytitle)
for n,C in enumerate(TSU.cumuIBD):
    C *= cxnorm
    TSP.cumuIBD[n] *= cxnorm

    l0, = ax.plot(xbins - dx0, C, label="%g s"%TSU.windows[n+1])
    ax.plot(xbins - dx0, TSP.cumuIBD[n], linestyle="--", color=l0.get_color())

l0, = ax.plot(xbins - 0.77, TSU.totalspect_ibd * cxnorm, label="total")
ax.plot(xbins - 0.77, TSP.totalspect_ibd * cxnorm, linestyle="--", color=l0.get_color())

ax.legend()
plt.ylim(0., None)
plt.xlim(0., 8)
fig.savefig("IBDs_time.pdf")

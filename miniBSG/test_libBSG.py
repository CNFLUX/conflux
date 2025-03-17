#!/bin/env python3

from ctypes import *
import numpy as np
import matplotlib.pyplot as plt
from conflux.BetaEngine import *

libBSG = cdll.LoadLibrary("libBSG.so")

#####################

libBSG.BSG_beta_spectrum.argtypes = [c_double, c_double, c_int, c_int, c_double, c_int, c_bool]
libBSG.BSG_beta_spectrum.restype = c_double

libBSG.BSG_FermiFunction.argtypes = [c_double, c_double, c_double]
libBSG.BSG_FermiFunction.restype = c_double

libBSG.BSG_AtomicScreeningCorrection.argtypes = [c_double, c_int]
libBSG.BSG_AtomicScreeningCorrection.restype = c_double

libBSG.BSG_RecoilCorrection.argtypes = [c_double, c_double, c_int]
libBSG.BSG_RecoilCorrection.restype = c_double

libBSG.BSG_QCorrection.argtypes = [c_double, c_double, c_int, c_int]
libBSG.BSG_QCorrection.restype = c_double

libBSG.BSG_L0Correction.argtypes = [c_double, c_int, c_double]
libBSG.BSG_L0Correction.restype = c_double

libBSG.BSG_RadiativeCorrection.argtypes = [c_double, c_double, c_int, c_double]
libBSG.BSG_RadiativeCorrection.restype = c_double

libBSG.shape_factor_unique_forbidden.argtypes = [c_double, c_int, c_double, c_int, c_double]
libBSG.shape_factor_unique_forbidden.restype = c_double

#####################


#####################
B = BetaBranch(55, 141, 0, 5., 5., 0., 1., 0.)
B.xbins = np.arange(0, B.Q, .05)

m_e = 0.51099895000
m_p = 938.27208816
m_n = m_p + 1.29333236
W0 = 1 + B.E0/m_e
neutron_R0 = 0.0025896*1.2
R = neutron_R0 * B.A**(1./3.)
M0 = (abs(B.Z)*m_p + (B.A-abs(B.Z))*m_n)/m_e
M2_F = 0
M2_GT = 1

# Conflux .py corrections
cfx = np.zeros(len(B.xbins))

cfx_ffn = []
cfx_fsz = []
cfx_scr = []
cfx_rgt = []
cfx_asc = []
cfx_rcx = []
cfx_cgt = []
cfx_sgt = []
cfx_sL0 = []

# BSG (wrapped) corrections
bsg = np.zeros(len(B.xbins))

bsg_ffn = []
bsg_fsz = []
bsg_rgt = []
bsg_rcx = []
bsg_asc = []
bsg_cgt = []

Ls = (2, 7, 11)
cfx_suf = {L: [] for L in Ls}
bsg_suf = {L: [] for L in Ls}

print("l = ", B.Parameters['l'])

np.seterr(divide='ignore', invalid='ignore', over='ignore')

for n,E in enumerate(B.xbins):
    W = 1 + E/m_e

    cfx[n] = B.BetaSpectrum(E)
    bsg[n] = libBSG.BSG_beta_spectrum(W, W0, B.A, B.Z, R, 0, False)

    #####
    bsg_ffn.append(libBSG.BSG_FermiFunction(W, B.Z, R))
    bsg_fsz.append(libBSG.BSG_L0Correction(W, B.Z, R))
    bsg_rcx.append(libBSG.BSG_RadiativeCorrection(W, W0, B.Z, R))
    bsg_asc.append(libBSG.BSG_AtomicScreeningCorrection(W,B.Z))

    bsg_cgt.append(libBSG.BSG_QCorrection(W, W0, B.Z, B.A))
    bsg_rgt.append(libBSG.BSG_RecoilCorrection(W, W0, B.A))

    #####
    cfx_ffn.append(fermi_function(W, **B.Parameters))
    cfx_fsz.append(finite_size_L0(W, **B.Parameters))

    cfx_rcx.append(radiative_correction(W, **B.Parameters))
    cfx_asc.append(atomic_screening(W, **B.Parameters))

    cfx_cgt.append(recoil_Coulomb_gamow_teller(W, **B.Parameters))
    cfx_rgt.append(recoil_gamow_teller(W, **B.Parameters))

    cfx_sgt.append(shape_factor_gamow_teller(W, **B.Parameters))
    cfx_sL0.append(finite_size_L0_simple(W, **B.Parameters))

    for L in Ls:
        bsg_suf[L].append(libBSG.shape_factor_unique_forbidden(W, L, W0, B.Z, R));
        cfx_suf[L].append(shape_factor_unique_forbidden(W, L, W0, B.Z, R));

#s_ratio = bsg*sum(cfx)/(cfx * sum(bsg))
s_ratio = bsg/cfx

if True:

    fig, ax = plt.subplots()
    ax.set(xlabel='E (MeV)', ylabel='beta decay spectrum')

    ax.plot(B.xbins, cfx, label="CONFLUX")
    ax.plot(B.xbins, bsg, label="BSG")
    ax.legend()

    fig.savefig("beta_spectrum.pdf")

if True:

    if True:
        plt.figure()
        fig, ax = plt.subplots()
        ax.set(xlabel='E (MeV)', ylabel='correction factor')

        ax.plot(B.xbins[1:], cfx_ffn[1:], label="CONFLUX Fermi Function")
        ax.plot(B.xbins[1:], bsg_ffn[1:], label="BSG Fermi Function (rescaled)", linestyle='dashed')

        ax.legend()
        fig.savefig("beta_Fermi.pdf")

    if True:
        plt.figure()
        fig, ax = plt.subplots()
        ax.set(xlabel='E (MeV)', ylabel='correction factor')

        ax.plot(B.xbins[1:], cfx_rcx[1:], label="CONFLUX radiative cxn")
        ax.plot(B.xbins[1:], bsg_rcx[1:], label="BSG radiative cxn", linestyle='dashed')

        ax.plot(B.xbins, cfx_sgt, label="CONFLUX GT shape", linestyle='dashed')

        ax.plot(B.xbins[1:], s_ratio[1:], label="BSG/CONFLUX")

        ax.legend()
        fig.savefig("beta_cxn.pdf")

    if True:
        plt.figure()
        fig, ax = plt.subplots()
        ax.set(xlabel='E (MeV)', ylabel='correction factor')

        for L in Ls:
            cfx_suf[L] = np.array(cfx_suf[L][1:])/sum(cfx_suf[L][1:])
            bsg_suf[L] = np.array(bsg_suf[L][1:])/sum(bsg_suf[L][1:])

            ax.plot(B.xbins[1:], cfx_suf[L], label="CONFLUX unique fbdn. L=%i"%L)
            ax.plot(B.xbins[1:], bsg_suf[L], label="BSG unique fbdn. L=%i"%L, linestyle='dashed')

        ax.legend()
        fig.savefig("beta_forbidden.pdf")

    if True:
        plt.figure()
        fig, ax = plt.subplots()
        ax.set(xlabel='E (MeV)', ylabel='correction factor')

        ax.plot(B.xbins, cfx_sL0, label="CONFLUX old L0 calc")
        ax.plot(B.xbins, cfx_fsz, label="CONFLUX L0 finite size")
        ax.plot(B.xbins[1:], bsg_fsz[1:], label="BSG L0 finite size", linestyle='dashed')

        ax.plot(B.xbins, cfx_asc, label="CONFLUX atomic screening")
        ax.plot(B.xbins, bsg_asc, label="BSG atomic screening", linestyle='dashed')

        # tiny
        # ax.plot(B.xbins, cfx_cgt, label="CONFLUX coulomb")
        #ax.plot(B.xbins, bsg_cgt, label="BSG coulomb", linestyle='dashed')

        #ax.plot(B.xbins, cfx_rgt, label="CONFLUX recoil")
        # ax.plot(B.xbins, bsg_rgt, label="BSG recoil", linestyle='dashed')

        ax.legend()
        fig.savefig("beta_matched.pdf")

import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import operator

# conflux modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine

e_fission = 3.2e-11 #joules
e_TNT = 4.184e9
r2 = 10000*2
fission_per_t = e_TNT/e_fission
test_flux = 15*fission_per_t/r2 # cm^-2

proton_per_cc = 6.81e22
ej309_dens = 0.962 # g_per_cc
detect_mass = 10e6 # g (e6 means ton)
proton_count = detect_mass/ej309_dens*proton_per_cc # per 10 ton

def ibd_xsection(enu):
    epos = enu - 1.15
    ppos = np.sqrt(epos**2 - 0.511*2)
    xsection = 0.0952*epos*ppos*1e-42 #Phys. Rev. D 60 053003 (Vogel model)
    xsection[np.isnan(xsection)] = 0
    return xsection

if __name__ == "__main__":
    xbins = np.arange(1.8, 15, 0.1)

    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB(defaultDB='JEFF')
    #U235.LoadCorrelation(defaultDB='ENDF')

    Pu239 = FissionIstp(94, 239)
    Pu239.LoadFissionDB(defaultDB='JEFF')
    # Pu239.LoadCorrelation()
    #U235.CalcCovariance(Ei=0)

    model = FissionModel()
    model.AddContribution(isotope=Pu239, Ei = 0.4, fraction=1, IFP=True)
    model.SaveToFile('FPY_239_JEFF_IFP_14MeV.csv')
    # model.AddContribution(isotope=Pu239, Ei = 0.4, fraction=1)
    # model.SaveToFile('FPY_239_JEFF_IFP.csv')

    #model.AddContribution(isotope=U233, Ei = 0, fraction=1)
    #model.AddContribution(isotope=Pu241, Ei = 0, fraction=0.0572)
    #model.AddIstp(39, 96, 1.0)

    # define the time windows of the calculation
    windows = np.logspace(-1, 6, 8)
    print('time:', windows)
    fig, ax = plt.subplots()
    spect_time = []
    print(ibd_xsection(xbins))

    # generate neutrino spectra in different time window after ignition
    for i in range(len(windows)-1):
        begin = windows[i]
        end = windows[i+1]

        sum_model= SumEngine(xbins = xbins)
        sum_model.AddModel(model, W=test_flux)

        betaSpectraDB = BetaEngine(sum_model.FPYlist.keys(), xbins=xbins)
        betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])
        sum_model.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False,  ifp_begin = begin,  ifp_end = end)
        spect = sum_model.spectrum
        print(np.argmax(spect), max(spect))
        spect_time.append(spect)

        ax.set(xlabel='E (MeV)', ylabel='neutrino/MeV', title='U-235 neutrino flux')
        ax.plot(sum_model.xbins, spect, label=str(begin)+' s - '+str(end)+' s')
    # ax.legend()
    # fig.savefig("239Pu_ENDF_jeff_0.4_MeV_time.png")

    # calcualte cumulative spectrum at time windows after the ignition
    fig, ax = plt.subplots()
    total_spect = np.zeros(len(xbins))
    i = 0
    y_time = []
    for spect in spect_time:
        print('spectrum intergral', sum(spect))
        total_spect+=spect
        ax.set(xlabel='E (MeV)', ylabel='neutrino/MeV', title='U-235 neutrino flux')
        ax.plot(xbins, total_spect, label='by '+ str(windows[i+1])+' s')
        i+=1
        y_time.append(sum(total_spect[xbins > 1.8]))
    fig.savefig('neutrino_overtime.png')


    fig, ax = plt.subplots()
    total_spect = np.zeros(len(xbins))
    i = 0
    y_time = []
    for spect in spect_time:
        print('spectrum intergral', sum(spect))
        total_spect+=spect
        ax.set(xlabel='E (MeV)', ylabel='IBD/MeV', title='U-235 neutrino flux')
        ax.plot(sum_model.xbins, total_spect*ibd_xsection(sum_model.xbins)*proton_count, label='by '+ str(windows[i+1])+' s')
        i+=1
        y_time.append(sum(total_spect[xbins > 1.8]))
    ax.legend()
    fig.savefig('IBD_overtime.png')
    # ax.legend()
    # fig.savefig("239Pu_ENDF_jeff_0.4_MeV_cumulative.png")
    print(windows[1:], y_time)
    fig, ax = plt.subplots()
    ax.set(xlabel='time (s)', ylabel='neutrino/MeV', title='evolution of neutrino flux')
    ax.plot(windows[1:], y_time)
    fig.savefig("test.png")

    spect_time = np.array(spect_time)
    ybins = windows
    print(ybins[1:])
    fig, ax = plt.subplots()
    pcmesh = ax.pcolormesh(xbins, ybins[1:], spect_time)
    ax.set_yscale('log')
    ax.set_ylim([1e-1, 1e7])
    fig.colorbar(pcmesh, ax=ax)
    fig.savefig("239Pu_ENDF_jeff_0.4_MeV_time_energy.png")

import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import operator

# conflux modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from scipy.interpolate import interp2d
import matplotlib.colors as colors

e_fission = 3.2e-11 #joules
e_TNT = 4.184e9
r2 = 10000**2
fission_per_t = e_TNT/e_fission
test_flux = 15*fission_per_t/r2 # cm^-2

proton_per_cc = 1e23
ej309_dens = 0.962 # g_per_cc
detect_mass = 1e3 # g (e6 means ton)
proton_count = detect_mass/ej309_dens*proton_per_cc # per 10 ton

det_eff = 0.5

def ibd_xsection(enu):
    epos = enu - 1.29
    epos[epos<0] = 0
    ppos = np.sqrt(epos**2-0.511**2)
    tau_n = 877.75
    fr = 1.7152
    xsection = 0.0952*epos*ppos*1e-42 #Phys. Rev. D 60 053003 (Vogel model) (m2)
    # xsection = 2*np.pi/0.511**5/(fr*tau_n)*epos*ppos #Phys. Rev. D 60 053003 (Vogel model) (m2)
    xsection[np.isnan(xsection)] = 0
    return xsection


if __name__ == "__main__":

    xbins = np.arange(0, 15, 0.1)

    U235 = FissionIstp(92, 235)
    U235.LoadFissionDB(DB='JEFF')
    #U235.LoadCorrelation(DB='ENDF')

    Pu239 = FissionIstp(94, 239)
    Pu239.LoadFissionDB(DB='JEFF')
    # Pu239.LoadCorrelation()
    #U235.CalcCovariance(Ei=0)

    model = FissionModel()
    model.AddContribution(isotope=U235, Ei = 0.4, fraction=1, IFP=True)
    model.SaveToFile('FPY_235_JEFF_IFP_14MeV.csv')
    # model.AddContribution(isotope=Pu239, Ei = 0.4, fraction=1)
    # model.SaveToFile('FPY_239_JEFF_IFP.csv')

    #model.AddContribution(isotope=U233, Ei = 0, fraction=1)
    #model.AddContribution(isotope=Pu241, Ei = 0, fraction=0.0572)
    #model.AddIstp(39, 96, 1.0)

    # define the time windows of the calculation
    windows = np.linspace(0, 100, 101)
    windows_log = np.logspace(-2, 4, 20)
    print('time:', windows)
    print('log time', windows_log)

    spect_time = []
    xsec = (ibd_xsection(xbins))
    plt.figure()
    plt.plot(xbins, xsec)
    plt.show()

    sum_model= SumEngine(xbins = xbins)
    sum_model.AddModel(model, W=test_flux)

    betaSpectraDB = BetaEngine(sum_model.FPYlist.keys(), xbins=xbins)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True, branchErange=[0.0, 20.0])
    sum_model.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False,  ifp_begin = 0,  ifp_end = 1e6)
    totalspect = sum_model.spectrum
    totalspect_ibd = sum_model.spectrum*ibd_xsection(xbins)*proton_count
    cumuspect = np.zeros(len(totalspect))
    fig, ax = plt.subplots()

    spect_sum_t = []
    spect_sum_t_thresh = []
    spect_sum_t_ibd = []

    # # generate neutrino spectra in different time window after ignition
    # for i in range(len(windows_log)-1):
    #     begin = windows_log[i]
    #     end = windows_log[i+1]
    #
    #     # sum_model.Clear()
    #     sum_model.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False,  ifp_begin = begin,  ifp_end = end)
    #     spect = sum_model.spectrum
    #     cumuspect += spect
    #     print(np.argmax(spect), max(spect))
    #     spect_time.append(spect)
    #     spect_sum_t.append(sum(cumuspect)/sum(totalspect)*100)
    #     spect_sum_t_thresh.append(sum(cumuspect[xbins > 1.8])/sum(totalspect[xbins > 1.8])*100)
    #     spect_sum_t_ibd.append(sum(cumuspect*ibd_xsection(xbins)*proton_count)/sum(totalspect_ibd)*100)
    #
    #     ax.set(xlabel='E (MeV)', ylabel='neutrino/MeV', title='U-235 neutrino flux')
    #     ax.plot(sum_model.xbins, spect, label=str(begin)+' s - '+str(end)+' s')
    #
    # ax.legend()
    # # plt.show()
    # fig.savefig("235U_ENDF_jeff_0.4_MeV_time.png")
    #
    # fig, ax = plt.subplots()
    # ax.set(xlabel='Time (s)', ylabel='Percent', title='U-235 neutrino flux over time')
    # ax.plot(windows_log[1:], spect_sum_t, label='total neutrino', color='red')
    # ax.plot(windows_log[1:], spect_sum_t_thresh, label='neutrino above 1.8 MeV', color='green')
    # ax.plot(windows_log[1:], spect_sum_t_ibd, label='IBDs', color='blue')
    # xval = np.array([1, 10, 100])
    # yval = np.interp(xval, windows_log[1:], spect_sum_t_ibd)
    # # Mark the data point on the curve
    # for j in range(len(xval)):
    #     ax.scatter(xval[j], yval[j], color='red')
    #     # Draw a line to the axes
    #     ax.axhline(y=yval[j], color='gray', linestyle='--', xmin=0, xmax=xval[j])
    #     ax.axvline(x=xval[j], color='gray', linestyle='--', ymin=0, ymax=yval[j])
    #
    # plt.grid(True, linestyle='--', alpha=0.7)
    # plt.legend()
    # plt.xscale('log')
    #
    # # plt.show()
    # fig.savefig("235U_ENDF_jeff_0.4_MeV_ratevstime_log.png")


    distance = np.logspace(0, 3, 20)
    tonage = np.logspace(-1, 2, 20)
    ibd_rate = []
    for iy in range(len(distance)):
        dista = distance[iy]
        # generate neutrino spectra in different time window after ignition
        ibd_rate_y = []
        for ix in range(len(tonage)):
            size = tonage[ix]

            r2 = (100*dista)**2
            test_flux = 15*fission_per_t/(4*np.pi*r2) # cm^-2

            detect_mass = 1e6*size # g (e6 means ton)
            proton_count = detect_mass/ej309_dens*proton_per_cc # per [size] ton

            # sum_model.Clear()
            sum_model= SumEngine(xbins = xbins)
            sum_model.AddModel(model, W=test_flux)
            sum_model.CalcReactorSpectrum(betaSpectraDB, branchErange=[0.0, 20.0], processMissing=False,  ifp_begin = 0,  ifp_end = 10)
            spect = sum_model.spectrum
            spect *= ibd_xsection(xbins)*proton_count*det_eff
            ibd_rate_y.append(sum(spect))
            print('tonage', size, 'distance', dista, sum(spect))

        ibd_rate.append(ibd_rate_y)

    # Interpolate the data

    ibd_rate = np.array(ibd_rate)
    fig, ax = plt.subplots()
    pcmesh = plt.pcolor(tonage, distance, ibd_rate,  norm=colors.LogNorm())

    # Add colorbar
    cbar = plt.colorbar()
    cbar.set_label('total IBD events')

    interp_func = interp2d(tonage, distance, ibd_rate, kind='linear')
    y_interp = np.logspace(0, 3, 100)
    x_interp = np.logspace(-1, 2, 100)

    z_interp = interp_func(x_interp, y_interp)
    # Save the interpolated data to a CSV file
    output_data = np.column_stack((x_interp, y_interp, z_interp))
    np.savetxt('interpolated_data1.csv', output_data, delimiter=',', header='X-axis,Y-axis,Interpolated_Value', comments='')

    # Plot the contour where interpolated values are equal to 1
    fig_contour = plt.contour(x_interp, y_interp, z_interp, levels=[1, 5, 10, 100], colors='red')
    x_locations = [20, 20, 20, 20]
    plt.clabel(fig_contour, inline=1, fontsize=8, inline_spacing=1)

    # Labeling axes
    plt.xlabel('Det mass (ton)')
    plt.ylabel('Distance (m)')
    ax.set_yscale('log')
    ax.set_xscale('log')

    # ax.set_ylim([1e-1, 1e7])
    # fig.colorbar(pcmesh, ax=ax)
    fig.savefig("235U_jeff_det_dist_0_10_eff.png")

    # # calcualte cumulative spectrum at time windows after the ignition
    # fig, ax = plt.subplots()
    # total_spect = np.zeros(len(xbins))
    # i = 0
    # y_time = []
    # for spect in spect_time:
    #     print('spectrum intergral', sum(spect))
    #     total_spect+=spect
    #     ax.set(xlabel='E (MeV)', ylabel='neutrino/MeV', title='U-235 neutrino flux')
    #     ax.plot(xbins, total_spect, label='by '+ str(windows[i+1])+' s')
    #     i+=1
    #     y_time.append(sum(total_spect[xbins > 1.8]))
    # fig.savefig('neutrino_overtime.png')
    #
    #
    # fig, ax = plt.subplots()
    # total_spect = np.zeros(len(xbins))
    # i = 0
    # y_time = []
    # for spect in spect_time:
    #     print('spectrum intergral', sum(spect))
    #     total_spect+=spect
    #     ax.set(xlabel='E (MeV)', ylabel='IBD/MeV', title='U-235 neutrino flux')
    #     ax.plot(sum_model.xbins, total_spect*ibd_xsection(sum_model.xbins)*proton_count, label='by '+ str(windows[i+1])+' s')
    #     i+=1
    #     y_time.append(sum(total_spect[xbins > 1.8]))
    # ax.legend()
    # fig.savefig('IBD_overtime.png')
    # # ax.legend()
    # # fig.savefig("239Pu_ENDF_jeff_0.4_MeV_cumulative.png")
    # print(windows[1:], y_time)
    # fig, ax = plt.subplots()
    # ax.set(xlabel='time (s)', ylabel='neutrino/MeV', title='evolution of neutrino flux')
    # ax.plot(windows[1:], y_time)
    # fig.savefig("test.png")

    # spect_time = np.array(spect_time)
    # ybins = windows
    # print(ybins[1:])
    # fig, ax = plt.subplots()
    # pcmesh = ax.pcolormesh(xbins, ybins[1:], spect_time)
    # ax.set_yscale('log')
    # ax.set_ylim([1e-1, 1e7])
    # fig.colorbar(pcmesh, ax=ax)
    # fig.savefig("235U_ENDF_jeff_0.4_MeV_time_energy.png")

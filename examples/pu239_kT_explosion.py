import numpy as np
import matplotlib.pyplot as plt

# conflux modules
from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp

import pandas as pd

# A simple burst fission model for fissions 100 ton TNT equivalent energy
e_fission = 3.2e-11 #joules
e_TNT = 4.184e9
r2 = 10000**2
fission_per_t = e_TNT/e_fission

det_eff = 0.5


# An IBD cross section empirical IBD cross section function  
# Phys. Rev. D 60 053003 (Vogel model) (m2)
def ibd_xsection(enu):
    epos = enu - 1.29
    epos[epos<0] = 0
    ppos = np.sqrt(epos**2-0.511**2)
    # tau_n = 877.75
    # fr = 1.7152
    xsection = 0.0952*epos*ppos*1e-42 
    # xsection = 2*np.pi/0.511**5/(fr*tau_n)*epos*ppos #Phys. Rev. D 60 053003 (Vogel model) (m2)
    xsection[np.isnan(xsection)] = 0
    return xsection


if __name__ == "__main__":

    xbins = np.arange(0, 10, 0.02)
    # define the time windows of the calculation
    windows = np.linspace(0, 100, 101)
    windows_log = np.logspace(-2, 3, 6)
    colorslist = ['red', 'blue', 'green', 'purple', 'orange', 'black']
    color_dict = dict(zip(windows_log, colorslist))
    print('time:', windows)
    print('log time', windows_log)
    
    # define lists of spectra in different time windows
    
    binwidth = xbins[1]-xbins[0]
    
    # Calculate beta spectra of all beta unstable isotopes
    betaSpectraDB = BetaEngine(xbins=xbins)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=True)

    # A simple burst fission model for fissions 100 ton TNT equivalent energy
    e_fission = 3.2e-11 #joules
    e_TNT = 4.184e9
    r2 = 10000**2
    fission_per_t = e_TNT/e_fission
    
    det_eff = 0.5
    
    distance = 50000   # cm2
    r2 = distance**2

    test_flux = 1 #1000*fission_per_t # cm^-2

    # Load the independent fission product yields from the JEFF database
    # Ei is set to be 0.4 MeV
    fission_istp_list = {}
    fission_istp_list["U233"] = FissionIstp(92, 233, 0.5, DB='ENDF', IFPY=True)
    fission_istp_list["U235"] = FissionIstp(92, 235, 0.5, DB='ENDF', IFPY=True)
    fission_istp_list["U238"] = FissionIstp(92, 238, 0.5, DB='ENDF', IFPY=True)
    fission_istp_list["Pu239"] = FissionIstp(94, 239, 0.5, DB='ENDF', IFPY=True)
    totalspect_list = {}
    totalibdspect_list = {}
    
    spect_time = {}
    spect_sum_t = {}
    spect_sum_t_thresh = {}
    spect_sum_t_ibd = {}
    average_e = {}
    neu_per_fission = {}
    neu_per_fission_per_sec = {}
    neu_per_fission_list = {}
    ibd_per_fission_list = {}
    
    for key, istp in fission_istp_list.items():
        istp.LoadFissionDB(DB='ENDF')
        istp.LoadCorrelation(DB='ENDF')
        istp.CalcBetaSpectra(betaSpectraDB)
        totalspect_list[key]=istp.spectrum*test_flux
        totalibdspect_list[key]=istp.spectrum*test_flux*ibd_xsection(xbins)

        cumuspect = np.zeros(len(xbins))
        
        spect_time[key] = []
        spect_sum_t[key] = []
        spect_sum_t_thresh[key] = []
        spect_sum_t_ibd[key] = []
        average_e[key] = []
        neu_per_fission[key] = []
        neu_per_fission_per_sec[key] = []
        neu_per_fission_list[key] = []
        ibd_per_fission_list[key] = []
        
        fig, ax = plt.subplots()


        # generate neutrino spectra in different time window after ignition
        for i in range(len(windows_log)-1):
            
            # calculate the beta spectra within different windows after the fission 
            begin = windows_log[i]
            end = windows_log[i+1]
            timelabel = str(begin)+' s - '+str(end)+' s'
            window_width = end-begin
            
            istp.CalcBetaSpectra(betaSpectraDB, processMissing=False,  ifp_begin = begin,  ifp_end = end)
            spect = istp.spectrum*test_flux
            
            average = sum(spect*xbins)/sum(spect)
            average_e[key].append(average)
            cumuspect += spect
            cumusum = sum(cumuspect)
            
            neu_per_fission[key].append(cumusum/test_flux*binwidth)
            neu_per_fission_per_sec[key].append(cumusum/test_flux*binwidth/window_width)
            spect_time[key].append(spect)
            spect_sum_t[key].append(cumusum/sum(totalspect_list[key])*100)  # total neutrino flux
            spect_sum_t_thresh[key].append(sum(cumuspect[xbins > 1.8])/sum(totalspect_list[key][xbins > 1.8])*100) #neutrinos above ibd threshold
            spect_sum_t_ibd[key].append(sum(cumuspect*ibd_xsection(xbins))/sum(totalibdspect_list[key])*100)    # measured IBD spectrum
        
            ibd_per_fission_list[key].append([timelabel, sum(cumuspect*ibd_xsection(xbins))/test_flux*binwidth])
            
            color = color_dict[begin]
            ax.plot(xbins, spect, label=timelabel, color=color)
        
        ax.set(xlabel='E (MeV)', ylabel='neutrino/MeV')
        # ax.plot(xbins, totalspect, label="total")
        ax.legend()
        ax.set_xlim(0, 10)
        ax.set_ylim(1e-6, 1)
        plt.yscale("log")
        fig.savefig(f"{key}_ENDF_jeff_0.4_MeV_time.pdf")
        
        neu_per_fission_list[key].append(["total", sum(totalspect_list[key])/test_flux*binwidth])
        ibd_per_fission_list[key].append(["total", sum(totalibdspect_list[key])/test_flux*binwidth])

        ibd_per_fission_df = pd.DataFrame(ibd_per_fission_list, columns=['time', key])
        ibd_per_fission_df.to_csv(f"{key}_ibd_per_fission.csv", index=False)
    
    fig, ax = plt.subplots()
    ax.set(xlabel='Time (s)', ylabel='neutrino IBD/fission/MeV')
    for key in totalibdspect_list:
        ax.plot(xbins, totalibdspect_list[key] , label=key)
    ax.set_xlim(1.8, 10)
    ax.legend()
    fig.savefig("ENDF_0.5_ibd_spect_per_fission.pdf")

    fig, ax = plt.subplots()
    ax.set(xlabel='Time (s)', ylabel='Average E (MeV)')
    ax.plot(windows_log[1:], average_e, label='Fast neutron', color='blue')
    plt.xscale('log')
    ax.legend()
    fig.savefig("239Pu_jeff_aveage_e.pdf")
    
    
    neu_per_fission_df = pd.DataFrame(neu_per_fission_list, columns=['time', 'Fast', '14 MeV'])
    neu_per_fission_df.to_csv("neu_per_fission.csv", index=False)
    
    fig, ax = plt.subplots()
    ax.set(xlabel='Time (s)', ylabel='neu/fission')
    ax.plot(windows_log[1:], neu_per_fission, label='Fast neutron', color='blue')
    ax.plot(windows_log[1:], neu_per_fission_14, label='14 MeV', color='red')
    ax.legend()
    plt.xscale('log')
    fig.savefig("239Pu_jeff_neu_per_fission.pdf")
    
    fig, ax = plt.subplots()
    ax.set(xlabel='Time (s)', ylabel='neu/fission/sec')
    ax.plot(windows_log[1:], neu_per_fission_per_sec, label='Fast neutron', color='blue')
    ax.plot(windows_log[1:], neu_per_fission_per_sec_14, label='14 MeV', color='red')
    plt.xscale('log')
    ax.legend()
    fig.savefig("239Pu_jeff_neu_per_sec.pdf")
    
    # Drawing the total neutrino flux with respect to time 
    fig, ax = plt.subplots()
    ax.set(xlabel='Time (s)', ylabel='Percent', title='Pu-239 neutrino flux over time (Fast)')
    ax.plot(windows_log[1:], spect_sum_t, label='total neutrino', color='red')
    ax.plot(windows_log[1:], spect_sum_t_thresh, label='neutrino above 1.8 MeV', color='green')
    ax.plot(windows_log[1:], spect_sum_t_ibd, label='IBDs', color='blue')
    xval = np.array([1, 10, 100])
    yval = np.interp(xval, windows_log[1:], spect_sum_t_ibd)
    # Mark the data point on the curve
    for j in range(len(xval)):
        ax.scatter(xval[j], yval[j], color='red')
        # Draw a line to the axes
        ax.axhline(y=yval[j], color='gray', linestyle='--', xmin=0, xmax=xval[j])
        ax.axvline(x=xval[j], color='gray', linestyle='--', ymin=0, ymax=yval[j])
    ax.set_ylim(0, 100)
    ax.legend()
    plt.xscale('log')
    fig.savefig("239Pu_ENDF_0.5_MeV_ratevstime_log.pdf")
    
    # Drawing the total neutrino flux with respect to time 
    fig, ax = plt.subplots()
    ax.set(xlabel='Time (s)', ylabel='Percent', title='Pu-239 neutrino flux over time (14 MeV)')
    ax.plot(windows_log[1:], spect_sum_t_14, label='total neutrino', color='red')
    ax.plot(windows_log[1:], spect_sum_t_thresh_14, label='neutrino above 1.8 MeV', color='green')
    ax.plot(windows_log[1:], spect_sum_t_ibd_14, label='IBDs', color='blue')
    xval = np.array([1, 10, 100])
    yval = np.interp(xval, windows_log[1:], spect_sum_t_ibd_14)
    # Mark the data point on the curve
    for j in range(len(xval)):
        ax.scatter(xval[j], yval[j], color='red')
        # Draw a line to the axes
        ax.axhline(y=yval[j], color='gray', linestyle='--', xmin=0, xmax=xval[j])
        ax.axvline(x=xval[j], color='gray', linestyle='--', ymin=0, ymax=yval[j])
    ax.set_ylim(0, 100)
    ax.legend()
    plt.xscale('log')
    fig.savefig("239Pu_ENDF_14_MeV_ratevstime_log.pdf")

    distance = np.logspace(5, 9, 20)
    tonage = np.logspace(3, 8, 20)
    ibd_rate = []
    
    test_flux = 15*fission_per_t # cm^-2
    ej309_dens = 0.962 # g_per_cc
    proton_per_cc = 3.16*1e23
    detect_mass = 1e3 # g (e6 means ton)
    proton_count = detect_mass/ej309_dens*proton_per_cc # per 10 ton
    
    Pu239.CalcBetaSpectra(betaSpectraDB, processMissing=False,  ifp_begin = 0,  ifp_end = 1000)
    spect = Pu239.spectrum

    for iy in range(len(distance)):
        dista = distance[iy]
        # generate neutrino spectra in different time window after ignition
        ibd_rate_y = []
        for ix in range(len(tonage)):
            size = tonage[ix]

            r2 = (100*dista)**2
            test_flux = 50000000*fission_per_t/(4*np.pi*r2) # cm^-2

            detect_mass = 1e6*size # g (e6 means ton)
            proton_count = detect_mass/ej309_dens*proton_per_cc # per [size] ton

            # sum_model.Clear()

            spect *= ibd_xsection(xbins)*proton_count*det_eff*test_flux
            ibd_rate_y.append(sum(spect))
            print('tonage', size, 'distance', dista, sum(spect))

        ibd_rate.append(ibd_rate_y)

    # Interpolate the data

    import matplotlib.colors as colors

    ibd_rate = np.array(ibd_rate)
    fig, ax = plt.subplots()
    pcmesh = plt.pcolor(tonage, distance, ibd_rate,  norm=colors.LogNorm())

    # Add colorbar
    cbar = plt.colorbar()
    cbar.set_label('total IBD events')

    from scipy.interpolate import RectBivariateSpline

    interp_func = RectBivariateSpline(tonage, distance, ibd_rate)
    y_interp = np.logspace(5, 9, 20)
    x_interp = np.logspace(3, 8, 20)

    z_interp = interp_func(x_interp, y_interp)  
    # Save the interpolated data to a CSV file
    output_data = np.column_stack((x_interp, y_interp, z_interp))
    np.savetxt('interpolated_data1.csv', output_data, delimiter=',', header='X-axis,Y-axis,Interpolated_Value', comments='')

    # Plot the contour where interpolated values are equal to 1
    fig_contour = plt.contour(x_interp, y_interp, z_interp, levels=[1, 10, 100, 1000], colors='red')
    x_locations = [20, 20, 20, 20]
    # plt.clabel(fig_contour, inline=1, fontsize=8, inline_spacing=1)

    # Labeling axes
    plt.xlabel('Det mass (ton)')
    plt.ylabel('Distance (m)')
    ax.set_yscale('log')
    ax.set_xscale('log')

    # ax.set_ylim([1e-1, 1e7])
    # fig.colorbar(pcmesh, ax=ax)
    fig.savefig("235U_jeff_det_dist_0_10_eff.pdf")

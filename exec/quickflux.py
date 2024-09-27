#!/usr/bin/env python3
import numpy as np
import json
import sys
import matplotlib.pyplot as plt
import csv

from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionModel, FissionIstp
from conflux.SumEngine import SumEngine
from conflux.ConversionEngine import ConversionEngine, BetaData


units = {'meV': 1e-9, 'eV': 1e-6, 'keV': 1e-3, 'MeV': 1}
JtoMeV = 6.242e12
rxunits = {"MW": 1e6*JtoMeV, "GW": 1e9*JtoMeV}

def read_json_file(json_file_path):
    """
    Reader of the jasonfile
    
    :param json_file_path: filename
    :type json_file_path: string
    :return: data
    :rtype: dict

    """
    
    try:
        # Read JSON data from the specified file path
        with open(json_file_path, "r") as file:
            data = json.load(file)
        return data
    except FileNotFoundError:
        print(f"Error: File '{json_file_path}' not found.")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in file '{json_file_path}'.")
        return None
    
def summation(json_file, xbins, nu_spec):
    summation = json_file['summation']
    qlow = summation['qlow']          # select the lower limit of the beta q values
    qhigh = summation['qhigh']        # select the upper limit of the beta q values
    missing = summation['missing']    # whether to process the missing beta isotopes
    IFP = ('time' in summation)
    ifp_begin = 0
    ifp_end = 0

    # Loading the beta spectrum data base and calculate beta spectra
    betaSpectraDB = BetaEngine(xbins=xbins)
    betaSpectraDB.CalcBetaSpectra(nu_spectrum=nu_spec, 
                                  branchErange=[qlow, qhigh])
    sum_model = SumEngine(betaSpectraDB)


    # Determine whether to calculate independent fission yield
    if IFP:
        "setting up reactor model dependent on time"
        ifptime = summation['time']

    #if 'fission_db' in json_file:
    
    # loading fission istps in the model
    fissionistps = summation['fissionistps']
    for fissile_istp in fissionistps:
        name = fissile_istp['name']
        Z = fissile_istp['Z']
        A = fissile_istp['A']
        Ei = fissile_istp['ei']
        istp = FissionIstp(Z, A, Ei, DB=fissile_istp['fissiondb'])
        if 'fissiondb' in fissile_istp:
            istp.LoadFissionDB(DB=fissile_istp['fissiondb'])
        if 'covariancedb' in fissile_istp:
            istp.LoadCovariance(DB=fissile_istp['covariancedb'])
        else:
            istp.LoadCovariance(DB=fissile_istp['fissiondb'])
        istp.CalcBetaSpectra(betaSpectraDB, 
                             processMissing=missing, 
                             ifp_begin = ifp_begin, 
                             ifp_end = ifp_end)

        count = fissile_istp['count']
        d_count = 0 if 'd_count' not in fissile_istp else fissile_istp["d_count"]
        sum_model.AddFissionIstp(istp, name, count, d_count)
    
    betaistps = summation['betaistps']
    for betaistps in fissionistps:
        name = betaistps['name']
        Z = betaistps['Z']
        A = betaistps['A']
    for istp in betaSpectraDB.istplist:
        betaIstp = betaSpectraDB.istplist[istp]
            
    sum_model.CalcReactorSpectrum()
        
    spectrum = sum_model.spectrum
    uncertainty = sum_model.uncertainty
    
    return spectrum, uncertainty

def conversion(json_file, xbins, nu_spec):
    conversion = json_file['conversion']
    if conversion:
        # Declare the conversion engine by adding beta data with corresponding 
        # FPY database
        convertmodel = ConversionEngine()
        fissionistps = conversion['fissionistps']
        for fissile_istp in fissionistps:
            name = fissile_istp['name']
            Z = fissile_istp['Z']
            A = fissile_istp['A']
            Ei = 0
            istp = FissionIstp(Z, A, Ei)
            istp.LoadFissionDB(DB="JEFF")

            beta_istp = BetaData(fissile_istp['conversiondb'])
            branch_slice = fissile_istp['slice']
            count = fissile_istp['count']
            d_count = fissile_istp['dcount']
            convertmodel.AddBetaData(beta_istp, istp, name, count, d_count=d_count)
            convertmodel.VBfitbeta(name, branch_slice)
        convert_spect, unc, cov = convertmodel.SummedSpectrum(xbins, 
                                                              nu_spectrum=nu_spec, 
                                                              cov_samp=5)
    
    return convert_spect, unc

def buildmodel(json_file):
    """
    Function to read items from the json file to setup a model and make run the
    calculation. 
    
    :param json_file: the input json data
    :type json_file: dict

    """
    # Setup spectrum range and binning
    xbins = np.arange(0, 20, 0.1)
    if 'spectrum' in json_file:
        spect = json_file['spectrum']
        xbins = np.arange(spect['binlow'], spect['binhigh'], spect['binwidth'])
    spectrum = np.zeros(len(xbins))
    uncertainty = np.zeros(len(xbins))

    # Setup beta spectrum calculation
    nu_spec = 1
    numass = 0
    if 'spectrum_type' in json_file:
        spectype = json_file['spectrum_type']
    # decide whether to calculate beta spectrum or neutrino spectrum
    if 'neutrino' in spectype:
        nu_spec = spectype['neutrino']
    # whether to include non-zero neutrino mass
    if 'numass' in spectype:
        numassarg = spectype['numass']
        numass = numassarg['value']*units[numassarg['unit']]

    # Reactor power setup for absolute neutrino flux calcualtion
    rxpower = 1
    if 'rxpower' in json_file:
        rxpower = json_file['rxpower']['value']*rxunits[json_file['rxpower']['unit']]

    # Calculate fissile isotopes with summation mode
    if 'summation' in json_file:
        sumspect, sumunc = summation(json_file, xbins, nu_spec)
        spectrum += sumspect
        uncertainty += sumunc

    # calculate the neutrino spectra with 
    if 'conversion' in json_file:
        convert_spect, unc = conversion(json_file, xbins, nu_spec)
        spectrum += convert_spect
        uncertainty += unc**2

    uncertainty = np.sqrt(uncertainty)
    
    # Stack the arrays together into a 2D array (transpose so each array is a column)
    data = np.column_stack((xbins, spectrum, uncertainty))
    
    # Save the array to a CSV file with column headers
    outputfile=json_file['output']
    np.savetxt(f'{outputfile}.csv', data, delimiter=',', header='e, spectrum, unc', comments='')
    
    fig, ax = plt.subplots()
    #ax.set_ylim([-1, 1])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', 
           ylabel='neutrino/decay/MeV', 
           title='U-233 neutrino flux')
    ax.errorbar(xbins, spectrum, yerr=uncertainty, label="w/ miss info")
    #ax.plot(sum2.xbins, summed_spect, label="w/o info")
    ax.legend()

    fig.savefig("233U_JEFF_quick.png")


if __name__ == "__main__":
    # Check if a JSON file path is provided as a command-line argument
    if len(sys.argv) < 2:
        print("Usage: python3 quickflux.py <json_file_path>")
        sys.exit(1)

    # Get the JSON file path from the command-line argument
    input_json = sys.argv[1]

    # Read JSON data from the specified file
    data = read_json_file(input_json)

    buildmodel(data)

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

def model_setup(json_file):

    # Setup spectrum range and binning
    xbins = np.arange(0, 20, 0.1)
    if 'spectrum' in json_file:
        spect = json_file['spectrum']
        xbins = np.arange(spect['binlow'], spect['binhigh'], spect['binwidth'])
    spectrum = np.zeros(len(xbins))

    # Setup beta spectrum calculation
    # decide whether to calculate beta spectrum or neutrino spectrum
    nu_spec = 1
    if 'neutrino' in json_file['beta_spec']:
        nu_spec = json_file['beta_spec']['neutrino']

    # whether to include non-zero neutrino mass
    numass = 0
    if 'numass' in json_file['beta_spec']:
        numassarg = json_file['beta_spec']['numass']
        numass = numassarg['value']*units[numassarg['unit']]

    # Reactor power setup for absolute neutrino flux calcualtion
    rxpower = 1
    if 'rxpower' in json_file:
        rxpower = json_file['rxpower']['value']*rxunits[json_file['rxpower']['unit']]

    # Calculate reactor composition with summation mode

    if 'sum_model' in json_file:
        summation = json_file['sum_model']
        qlow = summation['qlow']          # select the lower limit of the beta q values
        qhigh = summation['qhigh']        # select the upper limit of the beta q values
        missing = summation['missing']    # whether to process the missing beta isotopes
        IFP = ('time' in summation)
        ifp_begin = 0
        ifp_end = 0

        if IFP:
            "setting up reactor model dependent on time"
            ifptime = summation['time']

        # prepare a summation model
        fissionmodel = FissionModel()
        #if 'fission_db' in json_file:

        sum_composition = summation['composition']
        for fissile_istp in sum_composition:
            name = fissile_istp['name']
            Z = fissile_istp['Z']
            A = fissile_istp['A']
            Ei = fissile_istp['ei']
            istp = FissionIstp(fissile_istp['Z'], fissile_istp['A'])
            if 'fissiondb' in fissile_istp:
                istp.LoadFissionDB(defaultDB=fissile_istp['fissiondb'])
            if 'covariancedb' in fissile_istp:
                istp.LoadCovarience(defaultDB=fissile_istp['covariancedb'])

            fissionmodel.AddContribution(isotope=istp, Ei = Ei, fraction=1, IFP=IFP)

        sum_model= SumEngine(xbins = xbins)
        sum_model.AddModel(fissionmodel, W=1)
        betaSpectraDB = BetaEngine(sum_model.FPYlist.keys(), xbins=xbins)
        betaSpectraDB.CalcBetaSpectra(nu_spectrum=nu_spec, branchErange=[qlow, qhigh])
        for istp in betaSpectraDB.istplist:
            betaIstp = betaSpectraDB.istplist[istp]
            print(betaIstp.Q)
        sum_model.CalcReactorSpectrum(betaSpectraDB, branchErange=[qlow, qhigh], processMissing=missing,  ifp_begin = ifptime[0],  ifp_end = ifptime[1])
        spectrum += sum_model.spectrum

    conversion = json_file['convert_model']
    if conversion:
        # Declare the conversion engine by adding beta data with corresponding FPY
        # database
        convertmodel = ConversionEngine()
        conv_composition = conversion['composition']
        for fissile_istp in conv_composition:
            print(fissile_istp)
            name = fissile_istp['name']
            Z = fissile_istp['Z']
            A = fissile_istp['A']
            istp = FissionIstp(fissile_istp['Z'], fissile_istp['A'])
            istp.LoadFissionDB(defaultDB="JEFF")


            beta_istp = BetaData(fissile_istp['conversiondb'])
            branch_slice = fissile_istp['slice']
            fraction = fissile_istp['fraction']
            convertmodel.AddBetaData(beta_istp, istp, name, fraction)
            convertmodel.VBfitbeta(name, branch_slice)
        convert_spect, unc, cov = convertmodel.SummedSpectrum(xbins, nu_spectrum=True, cov_samp=5)
        spectrum += convert_spect

    print(spectrum)

    fig, ax = plt.subplots()
    #ax.set_ylim([-1, 1])
    #plt.yscale('log')
    ax.set(xlabel='E (MeV)', ylabel='neutrino/decay/MeV', title='U-235 neutrino flux')
    ax.plot(xbins, spectrum, label="w/ miss info")
    #ax.plot(sum2.xbins, summed_spect, label="w/o info")
    ax.legend()

    fig.savefig("233U_JEFF_quick.png")


if __name__ == "__main__":
    # Check if a JSON file path is provided as a command-line argument
    if len(sys.argv) < 2:
        print("Usage: python quickflux.py <json_file_path>")
        sys.exit(1)

    # Get the JSON file path from the command-line argument
    input_json = sys.argv[1]

    # Read JSON data from the specified file
    data = read_json_file(input_json)

    xbins = model_setup(data)
    print(xbins)

# Copyright 2025 Lawrence Livermore National Security, LLC. See the top-level NOTICE file for details.
# Author: Xianyi Zhang

# SPDX-License-Identifier: MIT

# fission product and spectra summation engine

"""Public modules"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm
import csv

"""CONFLUX modules."""
from conflux.config import CONFLUX_DB
from conflux.Basic import Spectrum

class NeutrinoData(Spectrum):
    """
    A class to save the beta spectrum data to be converted in the computer memory.

    ...

    Attributes
    ----------

    x : list
        A list of the x-axis bins
    y : list
        A list of the y-axis values (bin content)
    yerr : list
        A list of the errors on the y-axis (bin error)
    inputDB : str
        The filename where the spectrum data is saved
    spectrum : :class:`numpy.array`
        An array to store the full beta spectrum
    uncertainty : :class:`numpy.array`
        An array to store the full beta spectrum uncertainty


    Methods
    -------

    LoadConversionDB(self, inputDB, rel_err=True):
        Load the beta spectrum data
    """
    def __init__(self, inputDB, covarianceDB):
        self.inputDB = inputDB
        self.covarianceDB=covarianceDB
        
        self.LoadNeutrinoDB(self.inputDB, self.covarianceDB)

    def LoadNeutrinoDB(self, inputDB, covarianceDB):        
        df = pd.read_csv(inputDB,  comment='#')
        self.x = df["E"].to_numpy()
        self.y = df["Content"].to_numpy()

        Spectrum.__init__(self, self.x)

        self.spectrum = self.y
        self.cov_matrix = np.loadtxt(covarianceDB, delimiter=',')
        self.uncertainty = np.sqrt(self.cov_matrix.diagonal())

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    testnudata = []
    testnudata.append(NeutrinoData(f"{CONFLUX_DB}/neutrinoDB/U235_FLUX_DYB.csv", f"{CONFLUX_DB}/neutrinoDB/U235_FLUX_COV_DYB_ABS.csv"))
    testnudata.append(NeutrinoData(f"{CONFLUX_DB}/neutrinoDB/U235_FLUX_STEREO.csv", f"{CONFLUX_DB}/neutrinoDB/U235_FLUX_COV_STEREO.csv"))

    for i in range(len(testnudata)):
        plt.errorbar(testnudata[i].xbins, testnudata[i].spectrum, testnudata[i].uncertainty)
    # plt.ylim([0, 5e-43])
    plt.show()
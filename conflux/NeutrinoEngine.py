# Copyright 2025 Lawrence Livermore National Security, LLC. See the top-level NOTICE file for details.
# Author: Xianyi Zhang

# SPDX-License-Identifier: MIT

# fission product and spectra summation engine

"""Public modules"""
import numpy as np
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
        with open(os.path.expandvars(inputDB), newline='') as inputCSV:
            inputreader = csv.DictReader(inputCSV, delimiter=',', quotechar='|')